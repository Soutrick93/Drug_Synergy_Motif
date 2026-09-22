% =====================================================================
% run_all  -  Drug-synergy / network-motif pipeline: one entry point.
% ---------------------------------------------------------------------
% Place the public input files in <proj>/data (see data/README.md) and press
% Run. Every input is read from <proj>/data, every result is written to
% <proj>/results with a clear name, and every figure to <proj>/results/figures.
% Each step is a self-contained function file (helpers are local functions
% inside it), so there are no scattered one-line .m files to keep on the path.
%
% Steps:
%   step1_features : STRING + DGIdb + ALMANAC + SIGNOR  -> graph, pair features
%   step2_omics    : CellMiner expression + drug/cell PCA -> per-cell-line rows
%   step3_models   : network-only and molecular ablations, CV, statistics
%   step4_figures  : dataset, ablation and mechanism figures (PDF)
%   step5_baseline : propensity + fingerprint baseline comparison
%   step6_external : external validation on the O'Neil screen
%   step7_deep     : DeepSynergy-style MLP baseline on the same folds
% =====================================================================

clear; clc;

% -------------------- PROJECT LOCATION -------------------------------
% By default the project root is the folder that CONTAINS this code/ folder,
% so data/ and results/ are found automatically when you run this file from
% the repository. Override cfg.proj below if your data live elsewhere.
here     = fileparts(mfilename('fullpath'));   % .../drug-synergy-motifs/code
cfg.proj = fileparts(here);                    % .../drug-synergy-motifs
% cfg.proj = 'C:\path\to\drug-synergy-motifs'; % <-- uncomment to override
cfg.data    = fullfile(cfg.proj, 'data');
cfg.results = fullfile(cfg.proj, 'results');
cfg.figures = fullfile(cfg.results, 'figures');

% -------------------- INPUT FILE NAMES (in data/) --------------------
cfg.f_almanac      = 'almanac_synergy.csv';                 % ALMANAC/DrugComb ZIP
cfg.f_string_links = '9606.protein.links.v12.0.txt';        % STRING links
cfg.f_string_info  = '9606.protein.info.v12.0.txt';         % STRING protein->symbol
cfg.f_dgidb        = 'interactions.tsv';                     % DGIdb drug-target
cfg.f_signor       = 'all_data_15_09_26.tsv';               % SIGNOR causal network
cfg.f_expr         = 'nci60_expr.txt';                       % CellMiner NCI-60 expression

% -------------------- PARAMETERS -------------------------------------
cfg.string_thr    = 700;   % keep STRING edges with combined score >= this
cfg.pos_thr       = 10;    % ZIP >= this counts as synergistic in a cell line
cfg.neg_thr       = 0;     % ZIP <= this counts as non-synergistic (row level)
cfg.pos_minlines  = 3;     % pair-level positive: synergistic in >= this many lines
cfg.neg_maxlines  = 0;     % pair-level negative: synergistic in <= this many lines
cfg.dist_cap      = 12;    % cap on shortest target-target distance
cfg.k_drug        = 30;    % drug target-profile PCA dimensions
cfg.k_cell        = 30;    % cell expression PCA dimensions
cfg.K             = 5;     % cross-validation folds

% -------------------- RUN --------------------------------------------
if ~exist(cfg.results,'dir'), mkdir(cfg.results); end
if ~exist(cfg.figures,'dir'), mkdir(cfg.figures); end
addpath(fileparts(mfilename('fullpath')));      % so the step files are found

fprintf('==== drug-synergy pipeline ====\nproject: %s\n\n', cfg.proj);
step1_features(cfg);
step2_omics(cfg);
step3_models(cfg);
step4_figures(cfg);
step5_baseline(cfg);    % Figure 5 (needs no extra data; fingerprint bar optional)
step6_external(cfg);    % Figure 6 (needs data/oneil_synergy.csv; skips if absent)
step7_deep(cfg);        % deep-baseline comparison (needs fitcnet; optional)
fprintf('\n==== done. results in %s ====\n', cfg.results);
