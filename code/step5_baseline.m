function step5_baseline(cfg)
% STEP5_BASELINE  Compare the target-profile model against baselines under
% the same grouped CV, and draw Figure 5.
%   models: (1) drug-synergy propensity  (self-contained)
%           (2) target profile + cell    (this paper)
%           (3) chemical fingerprint + cell  (only if data/drug_fingerprints.csv
%               is present: name in column 1, fingerprint bits in the rest)
% Reads:  results/omics_features.mat  (+ optional data/drug_fingerprints.csv)
% Writes: results/baseline_results.mat, results/figures/Fig5_baseline.pdf
% Helper routines are local functions at the bottom of this file.

try, rng(1); catch, rand('seed',1); end
fprintf('=== step 5: baseline comparison ===\n');
BLUE=[0.165 0.471 0.839]; ORANGE=[0.922 0.408 0.204];
O = load(fullfile(cfg.results,'omics_features.mat'));
X=O.X; y=O.y; go=O.groups; rDA=O.rDA; rDB=O.rDB; rC=O.rC;
if isfield(O,'rZ'), rZ=O.rZ; else, error('omics_features.mat has no rZ; rerun step2_omics.'); end

% folds (same scheme as step3)
du=unique([rDA;rDB]); dmap=containers.Map(du,num2cell(mod(0:numel(du)-1,cfg.K)+1));
gA=cellfun(@(d)dmap(d),rDA); gB=cellfun(@(d)dmap(d),rDB);
cu=unique(rC); cmap=containers.Map(cu,num2cell(mod(0:numel(cu)-1,cfg.K)+1));
gC=cellfun(@(c)cmap(c),rC);
cellcols=find(go==2); ourcols=find(go<=2);

% optional fingerprint drug embedding, aligned to rows
[haveFP, fpA, fpB] = load_fp(cfg, rDA, rDB, cfg.k_drug);
mnames = {'drug propensity','target profile + cell'};
if haveFP, mnames{end+1} = 'fingerprint + cell'; end
nM = numel(mnames);

roc=cell(1,2); prc=cell(1,2);
for pm=1:2
    R=nan(cfg.K,nM); P=nan(cfg.K,nM);
    for k=1:cfg.K
        if pm==1, isTe=(gA==k)|(gB==k); isTr=(gA~=k)&(gB~=k);
        else,     isTe=(gC==k);          isTr=(gC~=k); end
        if sum(isTr)<50||sum(isTe)<10||numel(unique(y(isTr)))<2||numel(unique(y(isTe)))<2, continue; end
        % (1) propensity (train-only drug averages)
        Xp = propensity_features(rDA,rDB,rZ,isTr);
        R(k,1)=fit_eval(Xp(isTr,:),y(isTr),Xp(isTe,:),y(isTe),'roc');
        P(k,1)=fit_eval(Xp(isTr,:),y(isTr),Xp(isTe,:),y(isTe),'prc');
        % (2) target profile + cell (our model)
        R(k,2)=fit_eval(X(isTr,ourcols),y(isTr),X(isTe,ourcols),y(isTe),'roc');
        P(k,2)=fit_eval(X(isTr,ourcols),y(isTr),X(isTe,ourcols),y(isTe),'prc');
        % (3) fingerprint + cell
        if haveFP
            Xf=[fpA+fpB, abs(fpA-fpB), X(:,cellcols)];
            R(k,3)=fit_eval(Xf(isTr,:),y(isTr),Xf(isTe,:),y(isTe),'roc');
            P(k,3)=fit_eval(Xf(isTr,:),y(isTr),Xf(isTe,:),y(isTe),'prc');
        end
    end
    roc{pm}=R; prc{pm}=P;
end

% report
prot={'LEAVE-DRUG-OUT','LEAVE-CELL-LINE-OUT'};
for pm=1:2
    fprintf('\n%s\n', prot{pm});
    for s=1:nM
        [rm,rlo,rhi]=ci95(roc{pm}(:,s)); [pmv,plo,phi]=ci95(prc{pm}(:,s));
        fprintf('  %-22s AUROC %.3f [%.3f,%.3f]  AUPRC %.3f [%.3f,%.3f]\n', mnames{s},rm,rlo,rhi,pmv,plo,phi);
    end
end
save(fullfile(cfg.results,'baseline_results.mat'),'roc','prc','mnames','haveFP','-v7');

% -------- Figure 5: grouped bars, baseline vs model, per regime --------
base=mean(y==1);
mcols = {[0.66 0.70 0.74], BLUE, ORANGE};   % one colour per model
f=figure('Position',[80 80 860 340],'Color','w');
ax1=subplot(1,2,1); barpanel(ax1, roc, mnames, 0.5, 0.45, 1.0, '(a) AUROC', 'AUROC (95% CI)', mcols);
ax2=subplot(1,2,2); barpanel(ax2, prc, mnames, base, 0, 0.55, '(b) AUPRC', 'AUPRC (95% CI)', mcols);
savepdf(f, fullfile(cfg.figures,'Fig5_baseline'));
fprintf('\nwrote baseline_results.mat and Fig5_baseline.pdf\n');
if ~haveFP
    fprintf(['NOTE: no data/drug_fingerprints.csv found, so the chemical-fingerprint\n' ...
             '      bar is omitted. Add that file (see make_fingerprints.py) to include it.\n']);
end
end

% ===================== local: baselines ==============================
function Xp = propensity_features(rDA, rDB, rZ, trainmask)
    pm = containers.Map('KeyType','char','ValueType','double');
    cn = containers.Map('KeyType','char','ValueType','double');
    idx = find(trainmask(:))';
    for t = idx
        for dd = {rDA{t}, rDB{t}}
            d = dd{1};
            if isKey(pm,d), pm(d)=pm(d)+rZ(t); cn(d)=cn(d)+1; else, pm(d)=rZ(t); cn(d)=1; end
        end
    end
    ks = pm.keys; for i=1:numel(ks), pm(ks{i}) = pm(ks{i})/cn(ks{i}); end
    gmean = mean(rZ(trainmask));
    n = numel(rDA); Xp = zeros(n,2);
    for i = 1:n
        Xp(i,1) = getp(pm,rDA{i},gmean); Xp(i,2) = getp(pm,rDB{i},gmean);
    end
end
function v = getp(pm,d,g), if isKey(pm,d), v=pm(d); else, v=g; end, end

function [haveFP, fpA, fpB] = load_fp(cfg, rDA, rDB, kd)
    haveFP=false; fpA=[]; fpB=[];
    fpfile = fullfile(cfg.data,'drug_fingerprints.csv');
    if exist(fpfile,'file')~=2, return; end
    fid=fopen(fpfile,'r'); hdr=fgetl(fid); %#ok<NASGU>
    names={}; bits={};
    while true
        ln=fgetl(fid); if ~ischar(ln), break; end; if isempty(strtrim(ln)), continue; end
        c=strsplit(ln,','); names{end+1}=lower(strtrim(c{1}));
        bits{end+1}=cellfun(@str2double,c(2:end));
    end
    fclose(fid);
    FP=cell2mat(bits'); dmap=containers.Map(names,num2cell(1:numel(names)));
    kk=min(kd, min(size(FP))-1); Emb=pca_scores(FP,kk);
    n=numel(rDA); fpA=zeros(n,kk); fpB=zeros(n,kk);
    for i=1:n
        if isKey(dmap,rDA{i}), fpA(i,:)=Emb(dmap(rDA{i}),:); end
        if isKey(dmap,rDB{i}), fpB(i,:)=Emb(dmap(rDB{i}),:); end
    end
    haveFP=true;
end

function v = fit_eval(Xtr,ytr,Xte,yte,which)
    sc = train_predict(Xtr,ytr,Xte);
    if strcmp(which,'roc'), v=auroc(yte,sc); else, v=auprc(yte,sc); end
end

% ===================== local: models & metrics =======================
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
function [m,lo,hi]=ci95(x)
    x=x(~isnan(x)); nn=numel(x); m=mean(x); if nn<2, lo=m; hi=m; return; end
    tb=[12.706 4.303 3.182 2.776 2.571 2.447 2.365 2.306 2.262 2.228];
    df=nn-1; if df>=1&&df<=10, tc=tb(df); else, tc=1.96; end
    h=tc*std(x)/sqrt(nn); lo=m-h; hi=m+h;
end
function S=pca_scores(M,k)
    mu=mean(M,1); Mc=M-mu; [U,Sg,~]=svd(Mc,'econ'); sc=U*Sg; k=min(k,size(sc,2)); S=sc(:,1:k);
end

% ===================== local: figure =================================
function barpanel(ax, data, names, baseline, ylo, yhi, ttl, ylab, cols)
% grouped bars: two regime groups on x, one bar per model, 95% CI whiskers
    hold(ax,'on'); nM=numel(names); regimes={'leave-drug-out','leave-cell-line-out'};
    bw=0.8/nM; hleg=zeros(1,nM);
    for pm=1:2
        for s=1:nM
            [mm,lo,hi]=ci95(data{pm}(:,s));
            xc = pm + (s-(nM+1)/2)*bw; c=cols{s};
            hh=patch(ax, xc+[-bw/2 -bw/2 bw/2 bw/2], [0 mm mm 0], c, 'EdgeColor','w','LineWidth',0.5);
            if pm==1, hleg(s)=hh; end
            plot(ax,[xc xc],[lo hi],'-','Color',[0.28 0.28 0.27],'LineWidth',1.1);
            plot(ax,xc+[-0.05 0.05],[hi hi],'-','Color',[0.28 0.28 0.27],'LineWidth',1);
            plot(ax,xc+[-0.05 0.05],[lo lo],'-','Color',[0.28 0.28 0.27],'LineWidth',1);
        end
    end
    plot(ax,[0.4 2.6],[baseline baseline],'--','Color',[0.55 0.55 0.52],'LineWidth',1);
    set(ax,'XTick',1:2,'XTickLabel',regimes); xlim(ax,[0.4 2.6]); ylim(ax,[ylo yhi]);
    ylabel(ax,ylab); title(ax,ttl);
    legend(hleg, names, 'Location','southoutside','Orientation','horizontal','Box','off','FontSize',8);
    fig_style(ax);
end
function fig_style(ax)
    if nargin<1, ax=gca; end
    AX=[0.62 0.62 0.60]; TICKC=[0.20 0.20 0.19]; LAB=[0.05 0.05 0.05];
    set(ax,'FontName','Helvetica','FontSize',11,'LineWidth',0.75,'Box','off','TickDir','out', ...
        'TickLength',[0.012 0.012],'Layer','top','Color','w','XColor',AX,'YColor',AX);
    set(get(ax,'Title'),'FontSize',13,'FontWeight','bold','Color',LAB);
    set(get(ax,'XLabel'),'FontSize',12,'Color',LAB); set(get(ax,'YLabel'),'FontSize',12,'Color',LAB);
    try, ax.XAxis.TickLabelColor=TICKC; ax.YAxis.TickLabelColor=TICKC; catch, end
end
function savepdf(f, stem)
    drawnow; pdf=[stem '.pdf'];
    try, exportgraphics(f,pdf,'ContentType','vector');
    catch
        try, set(f,'PaperPositionMode','auto','Units','inches'); p=get(f,'Position');
             set(f,'PaperUnits','inches','PaperSize',[p(3) p(4)]); print(f,pdf,'-dpdf','-r300');
        catch, print(f,[stem '.png'],'-dpng','-r300'); end
    end
    try, savefig(f,[stem '.fig']); catch, end
end
