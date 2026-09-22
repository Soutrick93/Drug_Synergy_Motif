# The Marginal Value of Network Motifs in Anticancer Drug Synergy Prediction

Reproducible MATLAB pipeline for the study *"The Marginal Value of Network Motifs
in Anticancer Drug Synergy Prediction"* (Soutrick Das). The code tests whether
drug-pair network-motif features — in particular directed feed-forward loops —
predict anticancer drug synergy, using only public benchmark data and strict
grouped cross-validation.

The headline result is deliberately negative: motif counts do not beat a trivial
target-count baseline, because their apparent association with synergy is a
target-set-size (drug-promiscuity) confound. Molecular context (drug target
profiles plus cell-line expression) is what carries the signal, and the motifs
add only a marginal increment on top of it. The pipeline also benchmarks the
model against propensity, chemical-fingerprint, and DeepSynergy-style deep
baselines, and validates a frozen model on the independent O'Neil screen.

## What the pipeline produces

Running it end to end regenerates every quantitative result, table, and computed
figure in the paper:

- Pair-level network-motif ablation (Table 1)
- Per-cell-line molecular ablation with 95% confidence intervals (Table 2)
- Representation and model comparison — propensity, target profile, Morgan
  fingerprint, deep MLP (Table 3)
- Figures 1–5 as vector PDFs, plus the supplementary figures
- A plain-text statistics report (`results/stats_report.txt`)

## Repository layout

```
drug-synergy-motifs/
├── README.md              this file
├── LICENSE                MIT license
├── CITATION.cff           how to cite this software
├── code/                  the pipeline (MATLAB) + fingerprint helper (Python)
│   ├── run_all.m          one entry point: runs step1 → step7
│   ├── step1_features.m   build STRING+DGIdb+ALMANAC graph, SIGNOR net, pair features
│   ├── step2_omics.m      drug target-profile PCA + cell-expression PCA per row
│   ├── step3_models.m     grouped-CV ablations, statistics (Tables 1–2)
│   ├── step4_figures.m    dataset / ablation / mechanism figures (PDF)
│   ├── step5_baseline.m   propensity + chemical-fingerprint baselines (Table 3)
│   ├── step6_external.m   external validation on the O'Neil screen (Figure 5)
│   ├── step7_deep.m       DeepSynergy-style MLP baseline (Table 3)
│   └── make_fingerprints.py   SMILES → Morgan fingerprints (optional, needs RDKit)
├── data/                  place the public input files here (see data/README.md)
└── results/               all outputs are written here (created on first run)
    └── figures/           vector-PDF figures
```

## Requirements

- **MATLAB R2021a or later** with the **Statistics and Machine Learning Toolbox**
  (`fitcensemble` with the RUSBoost method for the tree ensembles; `fitcnet` for
  the deep baseline in `step7_deep.m`). No other toolboxes are required. If
  `fitcensemble` is unavailable the code falls back to a ridge-logistic scorer so
  the pipeline still runs, though the reported numbers assume RUSBoost.
- **Optional: Python 3 with RDKit** (`pip install rdkit`) only if you want to
  regenerate the chemical-fingerprint baseline from SMILES with
  `code/make_fingerprints.py`. A precomputed `data/drug_fingerprints.csv` can be
  used instead, and the fingerprint steps are skipped gracefully if it is absent.

## Data

**No datasets are redistributed in this repository.** Every input is publicly
available and must be downloaded from its original provider into the `data/`
folder. The datasets are governed by the licenses of their respective providers
and should be cited as in the paper. See **[`data/README.md`](data/README.md)**
for exact versions, filenames, expected columns, and references.

### Data sources

| File in `data/`                 | Dataset / version                        | Download |
|---------------------------------|------------------------------------------|----------|
| `almanac_synergy.csv`           | NCI-ALMANAC synergy (ZIP), via DrugComb  | https://drugcomb.org (study ALMANAC) |
| `9606.protein.links.v12.0.txt`  | STRING v12.0 (human) interactions        | https://string-db.org/cgi/download (taxon 9606) |
| `9606.protein.info.v12.0.txt`   | STRING v12.0 protein → gene symbol       | https://string-db.org/cgi/download (taxon 9606) |
| `interactions.tsv`              | DGIdb 5.0 drug–gene interactions         | https://www.dgidb.org/downloads |
| `all_data_15_09_26.tsv`         | SIGNOR causal network (human)            | https://signor.uniroma2.it/downloads.php |
| `nci60_expr.txt`                | CellMiner NCI-60 RNA-seq expression      | https://discover.nci.nih.gov/cellminer |
| `oneil_synergy.csv`             | O'Neil screen (optional, step 6)         | https://drugcomb.org (study ONEIL) |
| `drug_smiles.csv` → `drug_fingerprints.csv` | Morgan fingerprints (optional) | SMILES from https://pubchem.ncbi.nlm.nih.gov |

Primary references: Holbeck et al., *Cancer Res.* 2017 (ALMANAC); Zagidullin et
al., *NAR* 2019 and Zheng et al., *NAR* 2021 (DrugComb); Szklarczyk et al., *NAR*
2023 (STRING); Cannon et al., *NAR* 2024 (DGIdb); Licata et al., *NAR* 2020
(SIGNOR); Reinhold et al., *Cancer Res.* 2012/2019 (CellMiner); O'Neil et al.,
*Mol. Cancer Ther.* 2016 (O'Neil screen).

## How to run

1. Download the input files into `data/` (see `data/README.md`).
2. Open MATLAB, `cd` into the `code/` folder, and run:

   ```matlab
   run_all
   ```

   `run_all.m` sets the project root to the folder that contains `code/`, so
   `data/` and `results/` are found automatically. Override `cfg.proj` at the top
   of the file if your data live elsewhere.

3. Results appear in `results/` (`.mat` result files, `stats_report.txt`) and
   figures in `results/figures/`.

### Running a single step

Each step is a function that takes the same `cfg` struct. To run one step alone,
build `cfg` first, then call it. For example, to regenerate only the figures:

```matlab
cfg.proj    = fileparts(pwd);          % repo root, if you are in code/
cfg.data    = fullfile(cfg.proj,'data');
cfg.results = fullfile(cfg.proj,'results');
cfg.figures = fullfile(cfg.results,'figures');
cfg.f_expr  = 'nci60_expr.txt';
addpath(pwd);
step4_figures(cfg);
```

Steps 4–7 read the `.mat` files written by steps 1–3, so run `run_all` once (or
steps 1–3) before running a later step on its own.

## Key parameters

All parameters live at the top of `run_all.m` and are passed through `cfg`:

| Parameter        | Default | Meaning                                            |
|------------------|---------|----------------------------------------------------|
| `cfg.string_thr` | 700     | keep STRING edges with combined score ≥ this       |
| `cfg.pos_thr`    | 10      | ZIP ≥ this ⇒ synergistic (row level)               |
| `cfg.neg_thr`    | 0       | ZIP ≤ this ⇒ non-synergistic (row level)           |
| `cfg.pos_minlines` | 3     | pair positive if synergistic in ≥ this many lines  |
| `cfg.k_drug`     | 30      | drug target-profile PCA dimensions                 |
| `cfg.k_cell`     | 30      | cell-expression PCA dimensions                     |
| `cfg.K`          | 5       | cross-validation folds                             |

## Evaluation protocols

Two grouped cross-validation protocols are used throughout, both with five folds:

- **Leave-drug-out (LDO):** every held-out drug is absent from training.
- **Leave-cell-line-out (LCO):** every held-out cell line is absent from training.

Because synergy is rare (~5% positives), average precision (AUPRC) against the
base rate is the primary metric, reported alongside AUROC with 95% confidence
intervals. The fold-assignment scheme is identical across `step3`, `step5`, and
`step7`, so the baselines are compared on the same folds.

## Reproducibility note

The RUSBoost tree ensembles undersample the majority class stochastically. Each
step seeds the RNG (`rng(1)`), but the RNG state when a given model is fit depends
on how much code ran before it, so the same molecular model can differ by ≤0.002
AUROC/AUPRC between `step3` (Table 2) and `step5` (Table 3). This is expected and
does not affect any conclusion.

## Citing this software

If you use this code, please cite the paper and this repository (see
`CITATION.cff`).

## License

Released under the MIT License (see `LICENSE`). The input datasets are governed by
the licenses of their respective providers.
