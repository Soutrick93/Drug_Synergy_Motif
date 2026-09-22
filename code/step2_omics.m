function step2_omics(cfg)
% STEP2_OMICS  Build the per-cell-line feature matrix: drug target-profile
% PCA + cell expression PCA + the pair-level network features from step 1.
% Reads: results/graph.mat, results/pair_features.mat, data/<expr>
% Writes: results/omics_features.mat
% Helper routines are local functions at the bottom of this file.

fprintf('=== step 2: per-cell-line omics features ===\n');
G  = load(fullfile(cfg.results,'graph.mat'));
drug2tg = containers.Map(G.dkeys, G.dvals);
ndr = G.ndr; ndc = G.ndc; cl = G.cl; zip = G.zip; pairOK = G.pairOK;

% CellMiner expression
[cellNames, ~, expr] = read_cellminer(fullfile(cfg.data,cfg.f_expr));
fprintf('expression: %d genes x %d cell lines\n', size(expr,1), size(expr,2));

cl2col = containers.Map('KeyType','char','ValueType','double');
for j = 1:numel(cellNames), cl2col(normcell(cellNames{j})) = j; end

% drug target-profile embedding
drugs = {};
for i = 1:numel(ndr), if pairOK(i), drugs{end+1}=ndr{i}; drugs{end+1}=ndc{i}; end, end
drugs = unique(drugs); keep = false(1,numel(drugs));
for i = 1:numel(drugs), keep(i) = isKey(drug2tg,drugs{i}) && ~isempty(drug2tg(drugs{i})); end
drugs = drugs(keep);
tgU = {}; for i = 1:numel(drugs), tgU = [tgU, drug2tg(drugs{i})]; end
tgU = unique(tgU); tgIdx = containers.Map(tgU,num2cell(1:numel(tgU)));
DM = zeros(numel(drugs),numel(tgU));
for i = 1:numel(drugs), g = drug2tg(drugs{i}); for k = 1:numel(g), DM(i,tgIdx(g{k})) = 1; end, end
kd = min(cfg.k_drug, min(size(DM))-1);
drugEmb = pca_scores(DM, kd); drugRow = containers.Map(drugs,num2cell(1:numel(drugs)));
fprintf('drug embedding: %d drugs x %d dims\n', numel(drugs), kd);

% cell expression embedding (z-scored genes, PCA)
E = expr'; nanm = isnan(E);
if any(nanm(:))
    tmp = E; tmp(nanm) = 0; cnt = sum(~nanm,1); s = sum(tmp,1); cm = s./max(1,cnt);
    [~,ci] = find(nanm); E(nanm) = cm(ci);
end
mu = mean(E,1); sd = std(E,0,1); sd(sd==0) = 1; E = (E-mu)./sd;
kc = min(cfg.k_cell, size(E,1)-1); cellEmb = pca_scores(E, kc);
fprintf('cell embedding: %d cells x %d dims\n', size(cellEmb,1), kc);

% network group (pair features from step 1)
F = load(fullfile(cfg.results,'pair_features.mat'));
NX = F.X; netdim = size(NX,2);
netmap = containers.Map('KeyType','char','ValueType','double');
for p = 1:numel(F.pairA), sp = sort({F.pairA{p},F.pairB{p}}); netmap([sp{1} '|' sp{2}]) = p; end

% assemble rows
D = 2*kd + kc + netdim;
groups = [ones(1,2*kd), 2*ones(1,kc), 3*ones(1,netdim)];
maxN = numel(ndr); X = zeros(maxN,D); y = zeros(maxN,1); rZ = zeros(maxN,1);
rDA = cell(maxN,1); rDB = cell(maxN,1); rC = cell(maxN,1); nrow = 0;
for i = 1:numel(ndr)
    if ~pairOK(i) || isnan(zip(i)), continue; end
    a = ndr{i}; b = ndc{i};
    if ~isKey(drugRow,a) || ~isKey(drugRow,b), continue; end
    ck = normcell(cl{i}); if ~isKey(cl2col,ck), continue; end
    if zip(i)>=cfg.pos_thr, lab = 1; elseif zip(i)<=cfg.neg_thr, lab = 0; else, continue; end
    da = drugEmb(drugRow(a),:); db = drugEmb(drugRow(b),:); cvec = cellEmb(cl2col(ck),:);
    sp = sort({a,b}); nk = [sp{1} '|' sp{2}];
    if isKey(netmap,nk), nvec = NX(netmap(nk),:); else, nvec = zeros(1,netdim); end
    nrow = nrow+1; X(nrow,:) = [da+db, abs(da-db), cvec, nvec];
    y(nrow) = lab; rDA{nrow} = a; rDB{nrow} = b; rC{nrow} = cl{i}; rZ(nrow) = zip(i);
end
X = X(1:nrow,:); y = y(1:nrow); rDA = rDA(1:nrow); rDB = rDB(1:nrow); rC = rC(1:nrow); rZ = rZ(1:nrow);
fprintf('assembled rows: %d  (pos=%d, %.1f%% pos);  dims drug=%d cell=%d net=%d\n', ...
        nrow, sum(y==1), 100*mean(y==1), 2*kd, kc, netdim);
save(fullfile(cfg.results,'omics_features.mat'), 'X','y','groups','rDA','rDB','rC','rZ','-v7');
fprintf('wrote omics_features.mat\n\n');
end

% ===================== local helper functions ========================
function [cellNames, geneNames, expr] = read_cellminer(path)
    TAB = char(9); META = 6;
    fid = fopen(path,'r'); if fid<0, error('Cannot open %s', path); end
    cellNames = {}; ncell = 0;
    while true
        ln = fgetl(fid); if ~ischar(ln), fclose(fid); error('No header found'); end
        f = ssplit(ln,TAB);
        if ~isempty(f) && strncmpi(strtrim(f{1}),'Gene name',9)
            last = numel(f); while last>META && isempty(strtrim(f{last})), last = last-1; end
            cellNames = strtrim(f(META+1:last)); ncell = numel(cellNames); break;
        end
    end
    geneNames = {}; rows = {};
    while true
        ln = fgetl(fid); if ~ischar(ln), break; end
        if isempty(strtrim(ln)), continue; end
        f = ssplit(ln,TAB); if numel(f)<META+1, continue; end
        g = strtrim(f{1}); if isempty(g), continue; end
        v = nan(1,ncell);
        for c = 1:ncell, idx = META+c; if idx<=numel(f), v(c) = str2double(f{idx}); end, end
        geneNames{end+1,1} = g; rows{end+1,1} = v;
    end
    fclose(fid); expr = cell2mat(rows);
end

function parts = ssplit(ln, d)
    pos = [0, find(ln==d), numel(ln)+1]; parts = cell(1,numel(pos)-1);
    for k = 1:numel(pos)-1, parts{k} = ln(pos(k)+1:pos(k+1)-1); end
end

function c = normcell(s)
    c = upper(strtrim(s)); k = strfind(c,':'); if ~isempty(k), c = c(k(1)+1:end); end
    c = regexprep(c,'\(.*?\)',''); c = regexprep(c,'/ATCC',''); c = regexprep(c,'[^A-Z0-9]','');
end

function S = pca_scores(M, k)
    mu = mean(M,1); Mc = M - mu; [U,Sg,~] = svd(Mc,'econ');
    sc = U*Sg; k = min(k,size(sc,2)); S = sc(:,1:k);
end
