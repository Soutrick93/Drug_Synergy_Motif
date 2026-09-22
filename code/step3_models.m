function step3_models(cfg)
% STEP3_MODELS  Train and evaluate the network-only and molecular models
% under grouped cross-validation, with per-fold 95% CIs and a paired
% significance test for the network contribution.
% Reads:  results/pair_features.mat, results/omics_features.mat
% Writes: results/model_results.mat, results/stats_report.txt
% Helper routines are local functions at the bottom of this file.

try, rng(1); catch, rand('seed',1); end          % reproducible RUSBoost
fprintf('=== step 3: models, ablation, statistics ===\n');

% ---------------------------------------------------------------------
% A. Network-only ablation (pair level): leave-drug-out + random CV
% ---------------------------------------------------------------------
F = load(fullfile(cfg.results,'pair_features.mat'));
X = F.X; y = F.y; groups = F.groups; featnames = F.featnames; pairA = F.pairA; pairB = F.pairB;
lab = ~isnan(y); X = X(lab,:); y = y(lab); pairA = pairA(lab); pairB = pairB(lab);
drugs = unique([pairA, pairB]);
dfold = containers.Map('KeyType','char','ValueType','double');
for i = 1:numel(drugs), dfold(drugs{i}) = mod(i-1,cfg.K)+1; end
fA = cellfun(@(d)dfold(d),pairA); fB = cellfun(@(d)dfold(d),pairB);
sets = {find(groups==1), find(groups<=2), find(groups<=3)};
net_setname = {'B (target counts)','+STRING proximity','+SIGNOR FFL'};

net_roc = nan(cfg.K,3); net_prc = nan(cfg.K,3);
for k = 1:cfg.K
    isTe = (fA==k)|(fB==k); isTr = (fA~=k)&(fB~=k);
    for s = 1:3
        sc = train_predict(X(isTr,sets{s}),y(isTr),X(isTe,sets{s}));
        net_roc(k,s) = auroc(y(isTe),sc); net_prc(k,s) = auprc(y(isTe),sc);
    end
end
n = numel(y); perm = randperm(n); rfold = zeros(n,1); rfold(perm) = mod(0:n-1,cfg.K)+1;
net_rroc = nan(cfg.K,3); net_rprc = nan(cfg.K,3);
for k = 1:cfg.K
    isTe = (rfold==k); isTr = ~isTe;
    for s = 1:3
        sc = train_predict(X(isTr,sets{s}),y(isTr),X(isTe,sets{s}));
        net_rroc(k,s) = auroc(y(isTe),sc); net_rprc(k,s) = auprc(y(isTe),sc);
    end
end
net_imp = full_importance(X,y); net_base = mean(y==1);
fprintf('network ablation done (base AUPRC=%.3f)\n', net_base);

% ---------------------------------------------------------------------
% B. Molecular ablation (per cell line): leave-drug-out + leave-cell-out
% ---------------------------------------------------------------------
O = load(fullfile(cfg.results,'omics_features.mat'));
Xo = O.X; yo = O.y; go = O.groups; rDA = O.rDA; rDB = O.rDB; rC = O.rC;
osets = {find(go==1), find(go==2), find(go<=2), find(go<=3)};
om_setname = {'drug only','cell only','omics (drug+cell)','omics + network'};
du = unique([rDA;rDB]); dmap = containers.Map(du,num2cell(mod(0:numel(du)-1,cfg.K)+1));
gA = cellfun(@(d)dmap(d),rDA); gB = cellfun(@(d)dmap(d),rDB);
cu = unique(rC); cmap = containers.Map(cu,num2cell(mod(0:numel(cu)-1,cfg.K)+1));
gC = cellfun(@(c)cmap(c),rC); om_base = mean(yo==1);

om_roc = cell(1,2); om_prc = cell(1,2);   % {1}=LDO {2}=LCO ; each K x 4
for pm = 1:2
    R = nan(cfg.K,4); P = nan(cfg.K,4);
    for k = 1:cfg.K
        if pm==1, isTe = (gA==k)|(gB==k); isTr = (gA~=k)&(gB~=k);
        else,     isTe = (gC==k);          isTr = (gC~=k); end
        if sum(isTr)<50||sum(isTe)<10||numel(unique(yo(isTr)))<2||numel(unique(yo(isTe)))<2, continue; end
        for s = 1:4
            sc = train_predict(Xo(isTr,osets{s}),yo(isTr),Xo(isTe,osets{s}));
            R(k,s) = auroc(yo(isTe),sc); P(k,s) = auprc(yo(isTe),sc);
        end
    end
    om_roc{pm} = R; om_prc{pm} = P;
end
impo = full_importance(Xo,yo);
gimp = [sum(impo(go==1)) sum(impo(go==2)) sum(impo(go==3))]; gimp = gimp/sum(gimp);
fprintf('molecular ablation done (base AUPRC=%.3f)\n', om_base);

% ---------------------------------------------------------------------
% C. statistics report
% ---------------------------------------------------------------------
fid = fopen(fullfile(cfg.results,'stats_report.txt'),'w');
prot = {'LEAVE-DRUG-OUT','LEAVE-CELL-LINE-OUT'};
w2 = @(varargin) fprintf(fid, varargin{:});
w2('NETWORK-ONLY ABLATION (pair level, leave-drug-out) base AUPRC=%.3f\n', net_base);
for s = 1:3
    [rm,rlo,rhi]=ci95(net_roc(:,s)); [pmv,plo,phi]=ci95(net_prc(:,s));
    w2('  %-22s AUROC %.3f [%.3f,%.3f]  AUPRC %.3f [%.3f,%.3f]\n', net_setname{s},rm,rlo,rhi,pmv,plo,phi);
end
for pm = 1:2
    w2('\n=== %s ===  base AUPRC=%.3f\n', prot{pm}, om_base);
    for s = 1:4
        [rm,rlo,rhi]=ci95(om_roc{pm}(:,s)); [pmv,plo,phi]=ci95(om_prc{pm}(:,s));
        w2('  %-18s AUROC %.3f [%.3f,%.3f]  AUPRC %.3f [%.3f,%.3f]\n', om_setname{s},rm,rlo,rhi,pmv,plo,phi);
    end
    dR = om_roc{pm}(:,4)-om_roc{pm}(:,3); dP = om_prc{pm}(:,4)-om_prc{pm}(:,3);
    [dRm,dRlo,dRhi]=ci95(dR); [dPm,dPlo,dPhi]=ci95(dP);
    sR = tern(dRlo>0||dRhi<0,'significant','ns'); sP = tern(dPlo>0||dPhi<0,'significant','ns');
    w2('  network gain (omics+net - omics):\n');
    w2('    dAUROC=%.3f [%.3f,%.3f] %s ; dAUPRC=%.3f [%.3f,%.3f] %s\n', dRm,dRlo,dRhi,sR,dPm,dPlo,dPhi,sP);
    w2('    folds favouring +network: AUROC %d/%d, AUPRC %d/%d\n', sum(dR>0),sum(~isnan(dR)),sum(dP>0),sum(~isnan(dP)));
end
w2('\nfeature-group importance (molecular model): drug=%.2f cell=%.2f network=%.2f\n', gimp(1),gimp(2),gimp(3));
fclose(fid);
type_report(fullfile(cfg.results,'stats_report.txt'));

save(fullfile(cfg.results,'model_results.mat'), ...
     'net_roc','net_prc','net_rroc','net_rprc','net_setname','net_imp','net_base','featnames', ...
     'om_roc','om_prc','om_setname','om_base','gimp','-v7');
fprintf('wrote model_results.mat and stats_report.txt\n\n');
end

% ===================== local helper functions ========================
function score = train_predict(Xtr, ytr, Xte)
    mu = mean(Xtr,1); sd = std(Xtr,0,1); sd(sd==0)=1; Ztr = (Xtr-mu)./sd; Zte = (Xte-mu)./sd;
    if exist('fitcensemble','file')==2
        t = templateTree('MaxNumSplits',20);
        mdl = fitcensemble(Ztr,ytr,'Method','RUSBoost','NumLearningCycles',300,'Learners',t,'LearnRate',0.1);
        [~,sc] = predict(mdl,Zte); pos = find(mdl.ClassNames==1,1); score = sc(:,pos);
    else
        score = ridge_logistic(Ztr,ytr,Zte);
    end
end

function score = ridge_logistic(Ztr, ytr, Zte, lambda, iters)
    if nargin<4, lambda=1.0; end; if nargin<5, iters=50; end
    [n,d] = size(Ztr); Xd = [ones(n,1),Ztr]; yv = ytr(:); w = zeros(d+1,1);
    R = lambda*eye(d+1); R(1,1)=0;
    for it = 1:iters
        p = 1./(1+exp(-(Xd*w))); Wd = p.*(1-p)+1e-6;
        w = w - (Xd'*(Xd.*Wd)+R) \ (Xd'*(p-yv)+R*w);
    end
    score = 1./(1+exp(-[ones(size(Zte,1),1),Zte]*w));
end

function a = auroc(y, s)
    y = y(:); s = s(:); P = sum(y==1); N = sum(y==0);
    if P==0||N==0, a = NaN; return; end
    [ss,ord] = sort(s); yy = y(ord); n = numel(ss); ranks = zeros(n,1); i = 1;
    while i<=n
        j = i; while j<n && ss(j+1)==ss(i), j = j+1; end
        ranks(i:j) = (i+j)/2; i = j+1;
    end
    a = (sum(ranks(yy==1)) - P*(P+1)/2)/(P*N);
end

function ap = auprc(y, s)
    y = y(:); s = s(:); P = sum(y==1); if P==0, ap = NaN; return; end
    [~,ord] = sort(s,'descend'); yy = y(ord); tp = cumsum(yy==1); fp = cumsum(yy==0);
    prec = tp./(tp+fp); rec = tp/P; recprev = [0; rec(1:end-1)];
    ap = sum((rec-recprev).*prec);
end

function imp = full_importance(X, y)
    mu = mean(X,1); sd = std(X,0,1); sd(sd==0)=1; Z = (X-mu)./sd;
    if exist('fitcensemble','file')==2
        t = templateTree('MaxNumSplits',20);
        mdl = fitcensemble(Z,y,'Method','RUSBoost','NumLearningCycles',300,'Learners',t,'LearnRate',0.1);
        imp = predictorImportance(mdl);
    else
        n = size(Z,1); Xd = [ones(n,1),Z]; w = zeros(size(Xd,2),1); R = eye(size(Xd,2)); R(1,1)=0;
        for it = 1:50, p = 1./(1+exp(-(Xd*w))); Wd = p.*(1-p)+1e-6; w = w - (Xd'*(Xd.*Wd)+R)\(Xd'*(p-y)+R*w); end
        imp = abs(w(2:end))';
    end
end

function [m,lo,hi] = ci95(x)
    x = x(~isnan(x)); nn = numel(x); m = mean(x);
    if nn<2, lo = m; hi = m; return; end
    tb = [12.706 4.303 3.182 2.776 2.571 2.447 2.365 2.306 2.262 2.228];
    df = nn-1; if df>=1 && df<=10, tc = tb(df); else, tc = 1.96; end
    h = tc*std(x)/sqrt(nn); lo = m-h; hi = m+h;
end

function r = tern(c,a,b), if c, r = a; else, r = b; end, end

function type_report(p)
    fid = fopen(p,'r'); if fid<0, return; end
    while true, ln = fgetl(fid); if ~ischar(ln), break; end, disp(ln); end
    fclose(fid);
end
