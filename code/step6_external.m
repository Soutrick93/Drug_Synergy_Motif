function step6_external(cfg)
% STEP6_EXTERNAL  External validation on an independent screen (O'Neil).
% Trains the molecular (drug target-profile + cell expression) model on ALL
% of ALMANAC, then applies it -- frozen -- to the O'Neil pairs whose cell
% lines are shared with NCI-60, projecting O'Neil drugs and cells onto the
% SAME bases learned from ALMANAC. Draws Figure 6.
%
% Needs, in data/:  the O'Neil synergy table as  oneil_synergy.csv
%   columns (header, any order): drug_row, drug_col, cell_line_name, synergy_zip
%   (DrugComb study "ONEIL" export, or the O'Neil subset of DrugCombDB).
% Reads:  results/graph.mat (drug->target map + ALMANAC pairs), data/<expr>
% Writes: results/external_results.mat, results/figures/Fig6_external.pdf
% Helper routines are local functions at the bottom of this file.

try, rng(1); catch, rand('seed',1); end
fprintf('=== step 6: external validation (O''Neil) ===\n');
BLUE=[0.165 0.471 0.839]; ORANGE=[0.922 0.408 0.204]; MUTED=[0.537 0.529 0.506];
oneil = fullfile(cfg.data,'oneil_synergy.csv');
if exist(oneil,'file')~=2
    fprintf(['NOTE: data/oneil_synergy.csv not found. Download the O''Neil screen\n' ...
             '      (DrugComb study "ONEIL") with columns drug_row, drug_col,\n' ...
             '      cell_line_name, synergy_zip and place it in data/. Skipping.\n']);
    return;
end

% ---- ALMANAC drug->target map + expression basis --------------------
G = load(fullfile(cfg.results,'graph.mat'));
drug2tg = containers.Map(G.dkeys,G.dvals);
ndr=G.ndr; ndc=G.ndc; cl=G.cl; zip=G.zip; pairOK=G.pairOK;
[cellNames,~,expr] = read_cellminer(fullfile(cfg.data,cfg.f_expr));
cl2col = containers.Map('KeyType','char','ValueType','double');
for j=1:numel(cellNames), cl2col(normcell(cellNames{j})) = j; end

% drug target vocabulary + PCA basis, fit on ALMANAC drugs
adrugs = {}; for i=1:numel(ndr), if pairOK(i), adrugs{end+1}=ndr{i}; adrugs{end+1}=ndc{i}; end, end
adrugs = unique(adrugs); keep=false(1,numel(adrugs));
for i=1:numel(adrugs), keep(i)=isKey(drug2tg,adrugs{i})&&~isempty(drug2tg(adrugs{i})); end
adrugs = adrugs(keep);
tgU={}; for i=1:numel(adrugs), tgU=[tgU, drug2tg(adrugs{i})]; end
tgU=unique(tgU); tgIdx=containers.Map(tgU,num2cell(1:numel(tgU)));
DM = profile_matrix(adrugs, drug2tg, tgIdx);
kd = min(cfg.k_drug, min(size(DM))-1);
[drugEmb, muD, VD] = pca_fit(DM, kd);
drugRow = containers.Map(adrugs, num2cell(1:numel(adrugs)));

% cell expression basis (z-score genes, PCA), fit on all NCI-60 cells
E = expr'; nanm=isnan(E);
if any(nanm(:)), tmp=E; tmp(nanm)=0; cnt=sum(~nanm,1); s=sum(tmp,1); cm=s./max(1,cnt); [~,ci]=find(nanm); E(nanm)=cm(ci); end
muG=mean(E,1); sdG=std(E,0,1); sdG(sdG==0)=1; Ez=(E-muG)./sdG;
kc = min(cfg.k_cell, size(Ez,1)-1);
[cellEmb, muC, VC] = pca_fit(Ez, kc);

% ---- build ALMANAC training rows (drug + cell) ----------------------
[Xtr, ytr] = build_rows(ndr,ndc,cl,zip,pairOK, drugRow,drugEmb, cl2col,cellEmb, cfg);
fprintf('ALMANAC training rows: %d (%.1f%% pos)\n', numel(ytr), 100*mean(ytr==1));

% ---- read O'Neil and build external rows on shared cells ------------
S = read_named(oneil, ',');
odr = pick(S,{'drug_row','drugrow','drug1','druga','drug_1'});
odc = pick(S,{'drug_col','drugcol','drug2','drugb','drug_2'});
ocl = pick(S,{'cell_line_name','cellline','cell_line','cell','cell_line_id'});
[ozip, sctag] = pick(S,{'synergy_zip','zip','synergyzip','synergy_loewe','loewe','synergy_bliss','bliss','css','synergy','y'});
ozip = cellfun(@str2double, ozip);
% thresholds for the external metric (default to the training ZIP cutoffs)
pthr = getdef(cfg,'ext_pos_thr', cfg.pos_thr); nthr = getdef(cfg,'ext_neg_thr', cfg.neg_thr);
fprintf('external synergy column: %s   (positive if >= %.1f, negative if <= %.1f)\n', sctag, pthr, nthr);
odrn = cellfun(@canonical_drug, odr, 'UniformOutput',false);
odcn = cellfun(@canonical_drug, odc, 'UniformOutput',false);

Xte=[]; yte=[]; nmap=0; ncellmiss=0; ndrugmiss=0;
for i=1:numel(odrn)
    a=odrn{i}; b=odcn{i}; if isnan(ozip(i)), continue; end
    ck=normcell(ocl{i}); if ~isKey(cl2col,ck), ncellmiss=ncellmiss+1; continue; end
    da = project_drug(a, drug2tg, tgIdx, muD, VD, drugRow, drugEmb);
    db = project_drug(b, drug2tg, tgIdx, muD, VD, drugRow, drugEmb);
    if isempty(da)||isempty(db), ndrugmiss=ndrugmiss+1; continue; end
    if ozip(i)>=pthr, lab=1; elseif ozip(i)<=nthr, lab=0; else, continue; end
    cvec = cellEmb(cl2col(ck),:);
    Xte(end+1,:) = [da+db, abs(da-db), cvec]; yte(end+1,1)=lab; nmap=nmap+1;
end
fprintf('O''Neil rows usable: %d (cells not in NCI-60 skipped: %d; drugs unmapped: %d)\n', nmap, ncellmiss, ndrugmiss);
if nmap<20 || numel(unique(yte))<2
    fprintf('too few shared O''Neil rows for a stable estimate; stopping.\n'); return;
end

% ---- train on ALMANAC, score O'Neil ---------------------------------
sc = train_predict(Xtr, ytr, Xte);
AUC = auroc(yte,sc); AP = auprc(yte,sc); base=mean(yte==1);
fprintf('EXTERNAL (O''Neil): AUROC=%.3f  AUPRC=%.3f  (base rate=%.3f, n=%d, pos=%d)\n', AUC,AP,base,nmap,sum(yte==1));
save(fullfile(cfg.results,'external_results.mat'),'AUC','AP','base','nmap','yte','sc','-v7');

% ---- Figure 6: spider chart of metrics + confusion matrix -----------
% operating point: Youden-optimal threshold (maximises TPR - FPR)
P=sum(yte==1); N=sum(yte==0); ts=unique(sc); bestJ=-inf; tstar=median(sc);
for t=ts(:)'
    pred=sc>=t; tp=sum(pred&(yte==1)); fp=sum(pred&(yte==0));
    j=tp/max(1,P) - fp/max(1,N); if j>bestJ, bestJ=j; tstar=t; end
end
pred=sc>=tstar; TP=sum(pred&(yte==1)); FP=sum(pred&(yte==0)); FN=sum(~pred&(yte==1)); TN=sum(~pred&(yte==0));
prec_=TP/max(1,TP+FP); rec_=TP/max(1,TP+FN); spec=TN/max(1,TN+FP);
F1=2*prec_*rec_/max(1e-9,prec_+rec_); bacc=(rec_+spec)/2;
mcc=(TP*TN-FP*FN)/max(1e-9,sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)));
labels={'AUROC','AUPRC','Precision','Recall','F1','Specificity','Bal. acc.'};
vals  =[AUC, AP, prec_, rec_, F1, spec, bacc];

f=figure('Position',[80 80 900 380],'Color','w');
ax1=subplot(1,2,1); radar(ax1, labels, vals, BLUE);
title(ax1,'(a) External metrics (O''Neil)','FontWeight','bold','FontSize',13,'Color',[0.05 0.05 0.05]);
ax2=subplot(1,2,2); confmat(ax2, TP,FP,FN,TN);
title(ax2,sprintf('(b) Confusion matrix (MCC %.2f)',mcc),'FontWeight','bold','FontSize',13,'Color',[0.05 0.05 0.05]);
savepdf(f, fullfile(cfg.figures,'Fig6_external'));
fprintf('  operating point (Youden): precision=%.2f recall=%.2f F1=%.2f bal.acc=%.2f MCC=%.2f\n', prec_,rec_,F1,bacc,mcc);
fprintf('wrote external_results.mat and Fig6_external.pdf\n');
end

% ===================== local: feature construction ===================
function M = profile_matrix(drugs, drug2tg, tgIdx)
    M = zeros(numel(drugs), tgIdx.Count);
    for i=1:numel(drugs), g=drug2tg(drugs{i}); for k=1:numel(g), if isKey(tgIdx,g{k}), M(i,tgIdx(g{k}))=1; end, end, end
end

function [Emb, mu, V] = pca_fit(M, k)
    mu = mean(M,1); Mc = M - mu; [U,Sg,V] = svd(Mc,'econ');
    Emb = U*Sg; k=min(k,size(Emb,2)); Emb=Emb(:,1:k); V=V(:,1:k);
end

function e = project_drug(d, drug2tg, tgIdx, muD, VD, drugRow, drugEmb)
    e = [];
    if isKey(drugRow,d), e = drugEmb(drugRow(d),:); return; end   % seen in training
    if ~isKey(drug2tg,d), return; end
    v = zeros(1, tgIdx.Count); g = drug2tg(d); any1=false;
    for k=1:numel(g), if isKey(tgIdx,g{k}), v(tgIdx(g{k}))=1; any1=true; end, end
    if ~any1, return; end
    e = (v - muD) * VD;                                           % project onto frozen basis
end

function [X,y] = build_rows(ndr,ndc,cl,zip,pairOK, drugRow,drugEmb, cl2col,cellEmb, cfg)
    n=numel(ndr); X=zeros(n,2*size(drugEmb,2)+size(cellEmb,2)); y=zeros(n,1); m=0;
    for i=1:n
        if ~pairOK(i)||isnan(zip(i)), continue; end
        a=ndr{i}; b=ndc{i}; if ~isKey(drugRow,a)||~isKey(drugRow,b), continue; end
        ck=normcell(cl{i}); if ~isKey(cl2col,ck), continue; end
        if zip(i)>=cfg.pos_thr, lab=1; elseif zip(i)<=cfg.neg_thr, lab=0; else, continue; end
        da=drugEmb(drugRow(a),:); db=drugEmb(drugRow(b),:); cvec=cellEmb(cl2col(ck),:);
        m=m+1; X(m,:)=[da+db, abs(da-db), cvec]; y(m)=lab;
    end
    X=X(1:m,:); y=y(1:m);
end

function [col, name] = pick(S, cands)
    fn = fieldnames(S); col = {}; name = '';
    for i=1:numel(cands)
        for j=1:numel(fn)
            if strcmpi(fn{j}, cands{i}), col = S.(fn{j}); name = fn{j}; return; end
        end
    end
    error('column not found; looked for: %s', strjoin(cands,', '));
end

function v = getdef(s, f, d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% ===================== local: io, models, metrics, figure ============
function T = read_named(path, delim)
    fid=fopen(path,'r'); if fid<0, error('Cannot open %s', path); end
    hline=''; while true, ln=fgetl(fid); if ~ischar(ln), fclose(fid); error('no header'); end
        if ~isempty(strtrim(ln)), hline=ln; break; end, end
    heads=strsplit(hline,delim); nH=numel(heads); names=cell(1,nH);
    for i=1:nH, nm=regexprep(heads{i},'[^A-Za-z0-9_]',''); if isempty(nm), nm=sprintf('col%d',i); end
        if any(nm(1)=='0':'9'), nm=['x' nm]; end, names{i}=nm; end
    data=cell(0,nH); row=0;
    while true, ln=fgetl(fid); if ~ischar(ln), break; end; if isempty(strtrim(ln)), continue; end
        parts=strsplit(ln,delim,'CollapseDelimiters',false);
        if numel(parts)<nH, parts(numel(parts)+1:nH)={''}; end
        row=row+1; for i=1:nH, data{row,i}=strtrim(parts{i}); end
    end
    fclose(fid); T=struct(); for i=1:nH, if row==0, T.(names{i})={}; else, T.(names{i})=data(:,i); end, end
end

function s = normalize_drug(name)
    if isempty(name), s=''; return; end
    s=lower(strtrim(name)); s=regexprep(s,'\([^)]*\)',' '); s=regexprep(s,'\[[^\]]*\]',' '); s=strtrim(regexprep(s,'\s+',' '));
    salts={' hydrochloride',' dihydrochloride',' hydrobromide',' sulfate',' sulphate',' citrate',' mesylate', ...
           ' maleate',' tartrate',' ditartrate',' phosphate',' acetate',' succinate',' fumarate',' besylate', ...
           ' tosylate',' sodium',' potassium',' calcium',' chloride',' base',' anhydrous',' monohydrate',' dihydrate',' hcl',' hbr'};
    ch=true; while ch, ch=false; for k=1:numel(salts), L=numel(salts{k});
        if numel(s)>L && strcmp(s(end-L+1:end),salts{k}), s=strtrim(s(1:end-L)); ch=true; end, end, end
    s=strtrim(regexprep(s,'\s+',' '));
end

function s = canonical_drug(name)
    s=normalize_drug(name); persistent M
    if isempty(M), M=containers.Map('KeyType','char','ValueType','char');
        M('eloxatin')='oxaliplatin'; M('carboplatinum')='carboplatin'; M('cis-platin')='cisplatin';
        M('vepesid')='etoposide'; M('azacytidine')='azacitidine'; M('navelbine')='vinorelbine'; M('adm')='doxorubicin';
    end
    if isKey(M,s), s=M(s); end
end

function [cellNames, geneNames, expr] = read_cellminer(path)
    TAB=char(9); META=6; fid=fopen(path,'r'); if fid<0, error('Cannot open %s',path); end
    cellNames={}; ncell=0;
    while true, ln=fgetl(fid); if ~ischar(ln), fclose(fid); error('no header'); end
        f=ssplit(ln,TAB);
        if ~isempty(f)&&strncmpi(strtrim(f{1}),'Gene name',9)
            last=numel(f); while last>META&&isempty(strtrim(f{last})), last=last-1; end
            cellNames=strtrim(f(META+1:last)); ncell=numel(cellNames); break; end
    end
    geneNames={}; rows={};
    while true, ln=fgetl(fid); if ~ischar(ln), break; end; if isempty(strtrim(ln)), continue; end
        f=ssplit(ln,TAB); if numel(f)<META+1, continue; end; g=strtrim(f{1}); if isempty(g), continue; end
        v=nan(1,ncell); for c=1:ncell, idx=META+c; if idx<=numel(f), v(c)=str2double(f{idx}); end, end
        geneNames{end+1,1}=g; rows{end+1,1}=v; end
    fclose(fid); expr=cell2mat(rows);
end
function parts=ssplit(ln,d), pos=[0,find(ln==d),numel(ln)+1]; parts=cell(1,numel(pos)-1); for k=1:numel(pos)-1, parts{k}=ln(pos(k)+1:pos(k+1)-1); end, end
function c=normcell(s), c=upper(strtrim(s)); k=strfind(c,':'); if ~isempty(k), c=c(k(1)+1:end); end
    c=regexprep(c,'\(.*?\)',''); c=regexprep(c,'/ATCC',''); c=regexprep(c,'[^A-Z0-9]',''); end

function score = train_predict(Xtr, ytr, Xte)
    mu=mean(Xtr,1); sd=std(Xtr,0,1); sd(sd==0)=1; Ztr=(Xtr-mu)./sd; Zte=(Xte-mu)./sd;
    if exist('fitcensemble','file')==2
        t=templateTree('MaxNumSplits',20);
        mdl=fitcensemble(Ztr,ytr,'Method','RUSBoost','NumLearningCycles',300,'Learners',t,'LearnRate',0.1);
        [~,sc]=predict(mdl,Zte); pos=find(mdl.ClassNames==1,1); score=sc(:,pos);
    else
        [n,d]=size(Ztr); Xd=[ones(n,1),Ztr]; yv=ytr(:); w=zeros(d+1,1); R=eye(d+1); R(1,1)=0;
        for it=1:50, p=1./(1+exp(-(Xd*w))); Wd=p.*(1-p)+1e-6; w=w-(Xd'*(Xd.*Wd)+R)\(Xd'*(p-yv)+R*w); end
        score=1./(1+exp(-[ones(size(Zte,1),1),Zte]*w));
    end
end
function a=auroc(y,s)
    y=y(:); s=s(:); P=sum(y==1); N=sum(y==0); if P==0||N==0, a=NaN; return; end
    [ss,ord]=sort(s); yy=y(ord); n=numel(ss); ranks=zeros(n,1); i=1;
    while i<=n, j=i; while j<n&&ss(j+1)==ss(i), j=j+1; end, ranks(i:j)=(i+j)/2; i=j+1; end
    a=(sum(ranks(yy==1))-P*(P+1)/2)/(P*N);
end
function ap=auprc(y,s)
    y=y(:); s=s(:); P=sum(y==1); if P==0, ap=NaN; return; end
    [~,ord]=sort(s,'descend'); yy=y(ord); tp=cumsum(yy==1); fp=cumsum(yy==0);
    prec=tp./(tp+fp); rec=tp/P; recprev=[0; rec(1:end-1)]; ap=sum((rec-recprev).*prec);
end
function [fpr,tpr]=roc_curve(y,s)
    y=y(:); s=s(:); [~,ord]=sort(s,'descend'); yy=y(ord); P=sum(y==1); N=sum(y==0);
    tp=cumsum(yy==1); fp=cumsum(yy==0); tpr=[0; tp/max(1,P)]; fpr=[0; fp/max(1,N)];
end
function [rec,prec]=pr_curve(y,s)
    y=y(:); s=s(:); [~,ord]=sort(s,'descend'); yy=y(ord); P=sum(y==1);
    tp=cumsum(yy==1); fp=cumsum(yy==0); rec=tp/max(1,P); prec=tp./max(1,(tp+fp)); rec=[0;rec]; prec=[1;prec];
end
function radar(ax, labels, vals, col)
    axes(ax); cla(ax); hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
    n=numel(vals); ang=pi/2 - (0:n-1)/n*2*pi;
    th=linspace(0,2*pi,120);
    for r=0.2:0.2:1.0
        plot(ax, r*cos(th), r*sin(th), '-', 'Color',[0.86 0.86 0.84],'LineWidth',0.6);
        text(0.015, r, sprintf('%.1f',r),'FontSize',7,'Color',[0.62 0.62 0.60],'Parent',ax);
    end
    for i=1:n
        plot(ax,[0 cos(ang(i))],[0 sin(ang(i))],'-','Color',[0.86 0.86 0.84],'LineWidth',0.6);
        ha='center'; if cos(ang(i))>0.30, ha='left'; elseif cos(ang(i))<-0.30, ha='right'; end
        text(1.14*cos(ang(i)),1.14*sin(ang(i)),labels{i},'FontSize',9.5,'HorizontalAlignment',ha,'Color',[0.14 0.14 0.13],'Parent',ax);
    end
    px=vals(:)'.*cos(ang); py=vals(:)'.*sin(ang);
    patch('XData',[px px(1)],'YData',[py py(1)],'FaceColor',col,'FaceAlpha',0.20,'EdgeColor',col,'LineWidth',2,'Parent',ax);
    for i=1:n
        plot(ax,px(i),py(i),'o','MarkerFaceColor',col,'MarkerEdgeColor','w','MarkerSize',6);
        text(px(i)+0.04*cos(ang(i)),py(i)+0.04*sin(ang(i)),sprintf('%.2f',vals(i)),'FontSize',8,'Color',col,'FontWeight','bold','Parent',ax);
    end
    xlim(ax,[-1.4 1.4]); ylim(ax,[-1.3 1.35]);
end

function confmat(ax, TP,FP,FN,TN)
    axes(ax); M=[TP FN; FP TN];                     % rows: actual syn/non; cols: pred syn/non
    imagesc(ax, M); axis(ax,'square');
    cmap=[linspace(0.92,0.10,64)', linspace(0.95,0.34,64)', linspace(0.99,0.62,64)'];
    colormap(ax, cmap);
    lo=min(M(:)); hi=max(M(:));
    for r=1:2
        for c=1:2
            v=M(r,c); s=(v-lo)/max(1,(hi-lo)); tc=[0.1 0.1 0.1]; if s>0.6, tc=[1 1 1]; end
            text(c,r,num2str(v),'HorizontalAlignment','center','FontSize',16,'FontWeight','bold','Color',tc,'Parent',ax);
        end
    end
    set(ax,'XTick',[1 2],'XTickLabel',{'pred. syn','pred. non'}, ...
           'YTick',[1 2],'YTickLabel',{'actual syn','actual non'},'FontSize',10,'TickLength',[0 0]);
end

function fig_style(ax)
    if nargin<1, ax=gca; end; AX=[0.62 0.62 0.60]; TICKC=[0.20 0.20 0.19]; LAB=[0.05 0.05 0.05];
    set(ax,'FontName','Helvetica','FontSize',11,'LineWidth',0.75,'Box','off','TickDir','out','Layer','top','Color','w','XColor',AX,'YColor',AX);
    set(get(ax,'Title'),'FontSize',13,'FontWeight','bold','Color',LAB);
    set(get(ax,'XLabel'),'FontSize',12,'Color',LAB); set(get(ax,'YLabel'),'FontSize',12,'Color',LAB);
    try, ax.XAxis.TickLabelColor=TICKC; ax.YAxis.TickLabelColor=TICKC; catch, end
end
function savepdf(f, stem)
    drawnow; pdf=[stem '.pdf'];
    try, exportgraphics(f,pdf,'ContentType','vector');
    catch, try, set(f,'PaperPositionMode','auto','Units','inches'); p=get(f,'Position');
             set(f,'PaperUnits','inches','PaperSize',[p(3) p(4)]); print(f,pdf,'-dpdf','-r300');
        catch, print(f,[stem '.png'],'-dpng','-r300'); end, end
    try, savefig(f,[stem '.fig']); catch, end
end
