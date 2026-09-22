function step1_features(cfg)
% STEP1_FEATURES  Build the drug-target-PPI graph (STRING+DGIdb+ALMANAC),
% the directed signed SIGNOR network, and the pair-level feature table.
% Writes: results/graph.mat, results/signor.mat, results/pair_features.mat
% All helper routines are local functions at the bottom of this file.

TAB = char(9);
fprintf('=== step 1: build graph, SIGNOR network, and pair features ===\n');

% ---------------------------------------------------------------------
% 1. STRING protein -> gene symbol
% ---------------------------------------------------------------------
info = read_named(fullfile(cfg.data,cfg.f_string_info), TAB);
ensp2sym = containers.Map('KeyType','char','ValueType','char');
for i = 1:numel(info.string_protein_id)
    ensp2sym(info.string_protein_id{i}) = info.preferred_name{i};
end
fprintf('STRING info: %d proteins mapped\n', ensp2sym.Count);

% 2. STRING links -> symbol-level PPI adjacency
fid = fopen(fullfile(cfg.data,cfg.f_string_links),'r'); fgetl(fid);
C = textscan(fid,'%s %s %f'); fclose(fid);
p1 = C{1}; p2 = C{2}; sc = C{3};
keep = sc >= cfg.string_thr; p1 = p1(keep); p2 = p2(keep);
g1 = cell(numel(p1),1); g2 = cell(numel(p2),1); ok = false(numel(p1),1);
for i = 1:numel(p1)
    if isKey(ensp2sym,p1{i}) && isKey(ensp2sym,p2{i})
        g1{i}=ensp2sym(p1{i}); g2{i}=ensp2sym(p2{i}); ok(i)=true;
    end
end
g1 = g1(ok); g2 = g2(ok);
symbols = unique([g1; g2]); sym2idx = containers.Map(symbols,num2cell(1:numel(symbols)));
ei = cellfun(@(s)sym2idx(s),g1); ej = cellfun(@(s)sym2idx(s),g2);
m = numel(symbols); sel = ei~=ej; ei=ei(sel); ej=ej(sel);
A = sparse([ei;ej],[ej;ei],1,m,m); A = double(A>0);
fprintf('PPI graph: %d genes, %d undirected edges\n', m, nnz(A)/2);

% 3. DGIdb drug -> target genes
dg = read_named(fullfile(cfg.data,cfg.f_dgidb), TAB, true);
has_claim = isfield(dg,'drug_claim_name');
drug2tg = containers.Map('KeyType','char','ValueType','any');
for i = 1:numel(dg.drug_name)
    g = strtrim(dg.gene_name{i}); if isempty(g), continue; end
    keys = {}; d1 = normalize_drug(dg.drug_name{i}); if ~isempty(d1), keys{end+1}=d1; end
    if has_claim
        d2 = normalize_drug(dg.drug_claim_name{i});
        if ~isempty(d2) && (isempty(keys)||~strcmp(d2,keys{1})), keys{end+1}=d2; end
    end
    for kk = 1:numel(keys)
        d = keys{kk};
        if isKey(drug2tg,d), lst = drug2tg(d); else, lst = {}; end
        drug2tg(d) = unique([lst,{g}]);
    end
end
fprintf('DGIdb: %d normalised drug keys\n', drug2tg.Count);

% 4. ALMANAC synergy pairs
S = read_named(fullfile(cfg.data,cfg.f_almanac), ',');
dr = S.drug_row; dc = S.drug_col; cl = S.cell_line_name;
zip = cellfun(@str2double, S.synergy_zip);
ndr = cellfun(@canonical_drug, dr, 'UniformOutput', false);
ndc = cellfun(@canonical_drug, dc, 'UniformOutput', false);
npair = numel(dr);

pairOK = false(npair,1);
for i = 1:npair
    a = tg_in_graph(drug2tg,sym2idx,ndr{i});
    b = tg_in_graph(drug2tg,sym2idx,ndc{i});
    pairOK(i) = ~isempty(a) && ~isempty(b);
end
fprintf('ALMANAC rows: %d,  both drugs mapped: %d (%.0f%%)\n', ...
        npair, sum(pairOK), 100*sum(pairOK)/npair);

dkeys = drug2tg.keys(); dvals = drug2tg.values();
save(fullfile(cfg.results,'graph.mat'), 'symbols','A','dkeys','dvals', ...
     'ndr','ndc','cl','zip','pairOK','-v7');

% ---------------------------------------------------------------------
% 5. SIGNOR directed signed network
% ---------------------------------------------------------------------
R = read_named(fullfile(cfg.data,cfg.f_signor), TAB);
a = R.ENTITYA; b = R.ENTITYB;
keep = strcmp(R.TAX_ID,'9606') & strcmp(R.TYPEA,'protein') & strcmp(R.TYPEB,'protein');
for i = 1:numel(a)
    if keep(i) && (isempty(strtrim(a{i}))||isempty(strtrim(b{i}))), keep(i)=false; end
end
a = a(keep); b = b(keep); eff = R.EFFECT(keep);
sgn = zeros(numel(a),1);
for i = 1:numel(a)
    e = lower(strtrim(eff{i}));
    if strncmp(e,'up-regulates',12), sgn(i)=1; elseif strncmp(e,'down-regulates',14), sgn(i)=-1; end
end
dsym = unique([a;b]); sidx = containers.Map(dsym,num2cell(1:numel(dsym)));
si = cellfun(@(s)sidx(s),a); sj = cellfun(@(s)sidx(s),b);
nn = numel(dsym); sel = si~=sj; si=si(sel); sj=sj(sel); sgn=sgn(sel);
Ddir = double(sparse(si,sj,1,nn,nn)>0);
Dsig = spfun(@sign, sparse(si,sj,sgn,nn,nn));
sig2idx = containers.Map(dsym,num2cell(1:numel(dsym)));
fprintf('SIGNOR: %d nodes, %d directed edges (%d signed)\n', nn, nnz(Ddir), nnz(Dsig));
save(fullfile(cfg.results,'signor.mat'), 'dsym','Ddir','Dsig','-v7');

% ---------------------------------------------------------------------
% 6. pair-level features
% ---------------------------------------------------------------------
pkey = containers.Map('KeyType','char','ValueType','double');
pairA = {}; pairB = {}; zips = {};
for i = 1:numel(ndr)
    if ~pairOK(i) || isnan(zip(i)), continue; end
    da = ndr{i}; db = ndc{i}; if strcmp(da,db), continue; end
    sp = sort({da,db}); key = [sp{1} '|' sp{2}];
    if isKey(pkey,key), id = pkey(key);
    else, id = numel(pairA)+1; pkey(key)=id; pairA{id}=sp{1}; pairB{id}=sp{2}; zips{id}=[]; end
    zips{id}(end+1) = zip(i);
end
np = numel(pairA);
nSyn = zeros(np,1); nTested = zeros(np,1); medZip = zeros(np,1);
for p = 1:np
    nSyn(p)=sum(zips{p}>=cfg.pos_thr); nTested(p)=numel(zips{p}); medZip(p)=median(zips{p});
end
y = nan(np,1); y(nSyn>=cfg.pos_minlines)=1; y(nSyn<=cfg.neg_maxlines)=0;
fprintf('unordered mapped pairs: %d  (pos=%d neg=%d ambiguous=%d)\n', ...
        np, sum(y==1), sum(y==0), sum(isnan(y)));

% target index sets
drugs = unique([pairA, pairB]);
tgStr = containers.Map('KeyType','char','ValueType','any');
tgSig = containers.Map('KeyType','char','ValueType','any');
tgSym = containers.Map('KeyType','char','ValueType','any');
for i = 1:numel(drugs)
    d = drugs{i};
    if isKey(drug2tg,d), gs = drug2tg(d); else, gs = {}; end
    tgSym(d) = gs; is = []; ig = [];
    for k = 1:numel(gs)
        if isKey(sym2idx,gs{k}), is(end+1)=sym2idx(gs{k}); end
        if isKey(sig2idx,gs{k}), ig(end+1)=sig2idx(gs{k}); end
    end
    tgStr(d)=unique(is); tgSig(d)=unique(ig);
end

% STRING distances among target genes
allT = [];
for i = 1:numel(drugs), allT = [allT, tgStr(drugs{i})]; end
tgenes = unique(allT); tgpos = containers.Map(num2cell(tgenes),num2cell(1:numel(tgenes)));
Dtt = inf(numel(tgenes),numel(tgenes));
for i = 1:numel(tgenes), dvec = bfs_dist(A,tgenes(i)); Dtt(i,:) = dvec(tgenes)'; end

featnames = {'nTA','nTB','nShared','minDist','commonNbr','directEdge','unionNbr', ...
             'commonNbr_n','directEdge_n','coReg','convFFL','convCoh','convIncoh', ...
             'crossReg','reach2','coReg_n','convFFL_n','convIncoh_n','crossReg_n'};
groups = [1 1 1  2 2 2 2 2 2  3 3 3 3 3 3  3 3 3 3];
X = zeros(np, numel(featnames));
for p = 1:np
    da = pairA{p}; db = pairB{p};
    gaSym = tgSym(da); gbSym = tgSym(db);
    nTA = numel(gaSym); nTB = numel(gbSym); nShared = numel(intersect(gaSym,gbSym));
    As = tgStr(da); Bs = tgStr(db); denomStr = max(1, numel(As)*numel(Bs));
    if isempty(As) || isempty(Bs)
        minDist = cfg.dist_cap; commonNbr = 0; directEdge = 0; unionNbr = 0;
    else
        pa = cellfun(@(t)tgpos(t),num2cell(As)); pb = cellfun(@(t)tgpos(t),num2cell(Bs));
        d = min(min(Dtt(pa,pb))); if isinf(d), d = cfg.dist_cap; end
        minDist = d; NA = any(A(As,:),1); NB = any(A(Bs,:),1);
        commonNbr = full(sum(NA & NB)); directEdge = full(nnz(A(As,Bs))); unionNbr = full(sum(NA|NB));
    end
    [coReg,convFFL,convCoh,convIncoh,crossReg,reach2] = ffl_pair_features(tgSig(da),tgSig(db),Ddir,Dsig);
    denomSig = max(1, numel(tgSig(da))*numel(tgSig(db)));
    X(p,:) = [nTA,nTB,nShared, minDist,commonNbr,directEdge,unionNbr,commonNbr/denomStr,directEdge/denomStr, ...
              coReg,convFFL,convCoh,convIncoh,crossReg,reach2, coReg/denomSig,convFFL/denomSig,convIncoh/denomSig,crossReg/denomSig];
end
save(fullfile(cfg.results,'pair_features.mat'), 'X','y','featnames','groups', ...
     'pairA','pairB','nSyn','nTested','medZip','-v7');
fprintf('wrote graph.mat, signor.mat, pair_features.mat\n\n');
end

% ===================== local helper functions ========================
function T = read_named(path, delim, skipComments)
    if nargin < 3, skipComments = false; end
    fid = fopen(path,'r'); if fid<0, error('Cannot open file: %s', path); end
    hline = '';
    while true
        ln = fgetl(fid);
        if ~ischar(ln), fclose(fid); error('No header in: %s', path); end
        if isempty(strtrim(ln)), continue; end
        if skipComments && ~isempty(ln) && ln(1)=='#', continue; end
        hline = ln; break;
    end
    heads = local_split(hline, delim); n = numel(heads); names = cell(1,n);
    for i = 1:n
        nm = regexprep(heads{i}, '[^A-Za-z0-9_]', '');
        if isempty(nm), nm = sprintf('col%d',i); end
        if any(nm(1)=='0':'9'), nm = ['x' nm]; end
        names{i} = nm;
    end
    data = cell(0,n); row = 0;
    while true
        ln = fgetl(fid); if ~ischar(ln), break; end
        if isempty(strtrim(ln)), continue; end
        parts = local_split(ln, delim);
        if numel(parts) < n, parts(numel(parts)+1:n) = {''}; end
        row = row+1; for i = 1:n, data{row,i} = strtrim(parts{i}); end
    end
    fclose(fid); T = struct();
    for i = 1:n, if row==0, T.(names{i}) = {}; else, T.(names{i}) = data(:,i); end, end
end

function parts = local_split(ln, d)
    parts = {}; cur = ''; inq = false; i = 1; N = numel(ln);
    while i <= N
        c = ln(i);
        if inq
            if c=='"'
                if i<N && ln(i+1)=='"', cur(end+1)='"'; i=i+2; continue;
                else, inq=false; i=i+1; continue; end
            else, cur(end+1)=c; i=i+1; continue; end
        else
            if c=='"', inq=true; i=i+1; continue;
            elseif c==d, parts{end+1}=cur; cur=''; i=i+1; continue;
            else, cur(end+1)=c; i=i+1; continue; end
        end
    end
    parts{end+1} = cur;
end

function s = normalize_drug(name)
    if isempty(name), s = ''; return; end
    s = lower(strtrim(name));
    s = regexprep(s, '\([^)]*\)', ' '); s = regexprep(s, '\[[^\]]*\]', ' ');
    s = strtrim(regexprep(s, '\s+', ' '));
    salts = {' hydrochloride',' dihydrochloride',' hydrobromide',' hydroiodide', ...
             ' sulfate',' sulphate',' bisulfate',' hemisulfate',' citrate',' dicitrate', ...
             ' mesylate',' maleate',' malate',' tartrate',' ditartrate',' bitartrate', ...
             ' phosphate',' diphosphate',' acetate',' succinate',' fumarate',' besylate', ...
             ' tosylate',' pamoate',' lactate',' gluconate',' isethionate',' sodium', ...
             ' potassium',' calcium',' chloride',' bromide',' nitrate',' base',' anhydrous', ...
             ' monohydrate',' dihydrate',' trihydrate',' hcl',' hbr'};
    changed = true;
    while changed
        changed = false;
        for k = 1:numel(salts)
            L = numel(salts{k});
            if numel(s)>L && strcmp(s(end-L+1:end),salts{k}), s = strtrim(s(1:end-L)); changed = true; end
        end
    end
    s = strtrim(regexprep(s, '\s+', ' '));
end

function s = canonical_drug(name)
    s = normalize_drug(name);
    persistent M
    if isempty(M)
        M = containers.Map('KeyType','char','ValueType','char');
        M('eloxatin')='oxaliplatin'; M('carboplatinum')='carboplatin'; M('cis-platin')='cisplatin';
        M('emcyt')='estramustine'; M('trisenox')='arsenic trioxide'; M('zanosar')='streptozocin';
        M('vepesid j')='etoposide'; M('vepesid')='etoposide'; M('azacytidine')='azacitidine';
        M('azacytidine, 5-')='azacitidine'; M('5-fluoro-2''-deoxyuridine')='floxuridine';
        M('antibiotic ad 32')='valrubicin'; M('antibiotic ay 22989')='sirolimus';
        M('5-aminolevulinic acid')='aminolevulinic acid'; M('navelbine')='vinorelbine'; M('adm')='doxorubicin';
    end
    if isKey(M,s), s = M(s); end
end

function tg = tg_in_graph(drug2tg, sym2idx, d)
    tg = {}; if ~isKey(drug2tg,d), return; end
    cand = drug2tg(d);
    for k = 1:numel(cand), if isKey(sym2idx,cand{k}), tg{end+1} = cand{k}; end, end
end

function dist = bfs_dist(Aadj, src, maxd)
    n = size(Aadj,1); dist = inf(n,1); dist(src) = 0;
    frontier = false(n,1); frontier(src) = true;
    if nargin < 3, maxd = 25; end
    d = 0;
    while any(frontier) && d < maxd
        d = d+1; nb = any(Aadj(frontier,:),1)'; newf = nb & isinf(dist);
        dist(newf) = d; frontier = newf;
    end
end

function [coReg,convFFL,convCoh,convIncoh,crossReg,reach2] = ffl_pair_features(TA,TB,Ddir,Dsig)
    coReg=0; convFFL=0; convCoh=0; convIncoh=0; crossReg=0; reach2=0;
    if isempty(TA)||isempty(TB), return; end
    TA = TA(:)'; TB = TB(:)';
    outA = any(Ddir(TA,:),1); outB = any(Ddir(TB,:),1);
    coReg = full(sum(outA & outB)); crossReg = full(nnz(Ddir(TA,TB))+nnz(Ddir(TB,TA)));
    r2A = outA | any(Ddir(find(outA),:),1); r2B = outB | any(Ddir(find(outB),:),1);
    reach2 = full(sum(r2A & r2B));
    for x = TA
        for yv = TB
            if x==yv, continue; end
            if Ddir(x,yv)
                z = find(Ddir(x,:) & Ddir(yv,:)); z(z==x|z==yv) = [];
                for z1 = z, convFFL = convFFL+1; [convCoh,convIncoh] = tally(Dsig(x,yv),Dsig(yv,z1),Dsig(x,z1),convCoh,convIncoh); end
            end
            if Ddir(yv,x)
                z = find(Ddir(yv,:) & Ddir(x,:)); z(z==x|z==yv) = [];
                for z1 = z, convFFL = convFFL+1; [convCoh,convIncoh] = tally(Dsig(yv,x),Dsig(x,z1),Dsig(yv,z1),convCoh,convIncoh); end
            end
        end
    end
end

function [c,ic] = tally(sXY,sYZ,sXZ,c,ic)
    if sXY~=0 && sYZ~=0 && sXZ~=0
        if sign(sXY*sYZ)==sXZ, c = c+1; else, ic = ic+1; end
    end
end
