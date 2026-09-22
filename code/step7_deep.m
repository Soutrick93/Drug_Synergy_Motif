function step7_deep(cfg)
% STEP7_DEEP  Deep-baseline comparison. Trains a DeepSynergy-style multilayer
% perceptron on molecular features under the same two grouped protocols used
% throughout (leave-drug-out and leave-cell-line-out), so a published-style
% deep model is evaluated on the exact folds as the target-profile ensemble.
%
% Input features follow DeepSynergy's design (chemical representation plus
% cell-line expression): if data/drug_fingerprints.csv is present the drug is
% represented by Morgan-fingerprint principal components, otherwise by its
% target-profile principal components; cell lines are the expression principal
% components in both cases.
%
% Reads:  results/omics_features.mat (+ optional data/drug_fingerprints.csv)
% Writes: results/deep_results.mat and a line in results/stats_report.txt
%
% Requires MATLAB with fitcnet (Statistics and Machine Learning Toolbox).

    try, rng(1); catch, end
    fprintf('=== step 7: deep baseline ===\n');
    if exist('fitcnet','file') ~= 2
        error(['step7_deep needs fitcnet (Statistics and Machine Learning ' ...
               'Toolbox, MATLAB R2021a or later).']);
    end

    O = load(fullfile(cfg.results,'omics_features.mat'));
    X = O.X; y = O.y(:); go = O.groups; rDA = O.rDA; rDB = O.rDB; rC = O.rC;
    cellEmb = X(:, go==2);

    [haveFP, drugEmb, tag] = drug_representation(cfg, X, go, rDA, rDB);
    Z = [drugEmb, cellEmb];
    fprintf('deep baseline on %s + cell (%d features, %d rows)\n', tag, size(Z,2), size(Z,1));

    % grouped folds, identical scheme to steps 3 and 5
    du = unique([rDA;rDB]); dmap = containers.Map(du, num2cell(mod(0:numel(du)-1,cfg.K)+1));
    gA = cellfun(@(d)dmap(d), rDA); gB = cellfun(@(d)dmap(d), rDB);
    cu = unique(rC); cmap = containers.Map(cu, num2cell(mod(0:numel(cu)-1,cfg.K)+1));
    gC = cellfun(@(c)cmap(c), rC);

    prot = {'LDO','LCO'};
    roc = nan(cfg.K,2); prc = nan(cfg.K,2);
    for pm = 1:2
        for k = 1:cfg.K
            if pm==1, isTe = (gA==k)|(gB==k); isTr = (gA~=k)&(gB~=k);
            else,     isTe = (gC==k);          isTr = (gC~=k); end
            if sum(isTr)<50 || sum(isTe)<10 || numel(unique(y(isTr)))<2 || numel(unique(y(isTe)))<2
                continue;
            end
            sc = mlp_score(Z(isTr,:), y(isTr), Z(isTe,:));
            roc(k,pm) = auroc(y(isTe), sc);
            prc(k,pm) = auprc(y(isTe), sc);
        end
        [mr,lr,hr] = ci95(roc(:,pm)); [mp,lp,hp] = ci95(prc(:,pm));
        fprintf('  %s  AUROC %.3f [%.3f, %.3f]   AUPRC %.3f [%.3f, %.3f]\n', ...
                prot{pm}, mr, lr, hr, mp, lp, hp);
    end

    deep.tag = tag; deep.roc = roc; deep.prc = prc; deep.prot = prot;
    save(fullfile(cfg.results,'deep_results.mat'), '-struct', 'deep');

    fid = fopen(fullfile(cfg.results,'stats_report.txt'), 'a');
    if fid > 0
        fprintf(fid, '\nDeep baseline (%s + cell):\n', tag);
        for pm = 1:2
            [mr,lr,hr] = ci95(roc(:,pm)); [mp,lp,hp] = ci95(prc(:,pm));
            fprintf(fid, '  %s  AUROC %.3f [%.3f, %.3f]  AUPRC %.3f [%.3f, %.3f]\n', ...
                    prot{pm}, mr, lr, hr, mp, lp, hp);
        end
        fclose(fid);
    end
    fprintf('wrote deep_results.mat\n');
end

% ===================== local helper functions ========================
function sc = mlp_score(Xtr, ytr, Xte)
% DeepSynergy-style feed-forward network: two hidden layers, ReLU, input
% standardisation, trained on the raw class balance and scored by the
% positive-class posterior.
    net = fitcnet(Xtr, ytr, 'Standardize', true, ...
        'LayerSizes', [512 256], 'Activations', 'relu', ...
        'Lambda', 1e-4, 'IterationLimit', 300, 'Verbose', 0);
    [~, post] = predict(net, Xte);
    pc = find(net.ClassNames == 1, 1);
    if isempty(pc), pc = size(post,2); end
    sc = post(:, pc);
end

function [haveFP, emb, tag] = drug_representation(cfg, X, go, rDA, rDB)
% Morgan-fingerprint drug embedding if data/drug_fingerprints.csv is present,
% otherwise the target-profile drug embedding already in the feature matrix.
    haveFP = false; emb = X(:, go==1); tag = 'target profile';
    fpfile = fullfile(cfg.data, 'drug_fingerprints.csv');
    if ~exist(fpfile,'file'), return; end
    T = readtable(fpfile, 'ReadVariableNames', true);
    names = lower(strtrim(string(T{:,1})));
    bits = table2array(T(:, 2:end));
    fp = containers.Map('KeyType','char','ValueType','any');
    for i = 1:numel(names), fp(char(names(i))) = bits(i,:); end
    nb = size(bits,2);
    getfp = @(d) local_getfp(fp, lower(d), nb);
    FA = cell2mat(cellfun(@(d) getfp(d), rDA, 'UniformOutput', false));
    FB = cell2mat(cellfun(@(d) getfp(d), rDB, 'UniformOutput', false));
    pair = [FA + FB, abs(FA - FB)];          % order-invariant pair encoding
    k = min(cfg.k_drug, min(size(pair)) - 1);
    mu = mean(pair,1); Pc = pair - mu; [U,S,~] = svd(Pc, 'econ');
    emb = U(:,1:k) * S(1:k,1:k);
    haveFP = true; tag = 'Morgan fingerprint';
end

function v = local_getfp(fp, d, nb)
    if isKey(fp, d), v = fp(d); else, v = zeros(1, nb); end
end

function a = auroc(y, s)
    y = y(:) == 1; s = s(:); [~, o] = sort(s, 'descend'); y = y(o);
    P = sum(y); N = sum(~y); if P==0 || N==0, a = 0.5; return; end
    tp = cumsum(y); fp = cumsum(~y);
    a = trapz([0; fp/N], [0; tp/P]);
end

function ap = auprc(y, s)
    y = y(:) == 1; s = s(:); [~, o] = sort(s, 'descend'); y = y(o);
    tp = cumsum(y); fp = cumsum(~y); P = sum(y); if P==0, ap = 0; return; end
    prec = tp ./ max(1, tp + fp); rec = tp / P;
    ap = sum(diff([0; rec]) .* prec);
end

function [m, lo, hi] = ci95(x)
    x = x(~isnan(x)); n = numel(x); m = mean(x);
    if n < 2, lo = m; hi = m; return; end
    tb = [12.706 4.303 3.182 2.776 2.571 2.447 2.365 2.306 2.262 2.228];
    df = n - 1; if df >= 1 && df <= 10, tc = tb(df); else, tc = 1.96; end
    h = tc * std(x) / sqrt(n); lo = m - h; hi = m + h;
end
