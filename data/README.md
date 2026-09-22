# Input data

No data are redistributed here. Every input is public and must be downloaded from
the source below and placed in this `data/` folder with the exact filename shown.
The datasets are governed by the licenses of their respective providers; cite them
as in the paper.

The pipeline maps 89 of 103 ALMANAC drugs to network targets and 59 of 60 cell
lines to CellMiner; exact counts depend on the database releases you download.

## Required files

### `almanac_synergy.csv` — NCI-ALMANAC synergy (ZIP), via DrugComb
- Source: DrugComb portal, https://drugcomb.org (study **ALMANAC**), or the
  NCI-ALMANAC release harmonised through DrugComb.
- Content: drug-pair / cell-line synergy measurements with a ZIP synergy score.
- Expected columns (header, any order): a drug-1 name, a drug-2 name, a cell-line
  name, and a ZIP synergy value. `step1_features.m` reads these by name; open the
  file once to confirm the column headers match what step 1 expects and adjust the
  header names there if your export differs.
- Reference: Holbeck et al., *Cancer Res.* 2017; Zagidullin et al., *NAR* 2019;
  Zheng et al., *NAR* 2021.

### `9606.protein.links.v12.0.txt` and `9606.protein.info.v12.0.txt` — STRING v12.0
- Source: https://string-db.org/cgi/download (organism *Homo sapiens*, taxon
  9606), version 12.0.
  - `9606.protein.links.v12.0.txt` — protein–protein links with combined score.
  - `9606.protein.info.v12.0.txt` — protein ID → preferred gene symbol.
- The pipeline keeps edges with combined score ≥ `cfg.string_thr` (default 700).
- Reference: Szklarczyk et al., *NAR* 2023.

### `interactions.tsv` — DGIdb drug–gene interactions
- Source: https://www.dgidb.org/downloads (interactions TSV), DGIdb 5.0.
- Content: drug → gene (target) interactions; used to build each drug's
  protein-target set.
- Reference: Cannon et al., *NAR* 2024.

### `all_data_15_09_26.tsv` — SIGNOR causal network
- Source: https://signor.uniroma2.it/downloads.php (full human dump, TSV).
- Content: directed, signed causal edges (activation / inhibition) used for the
  feed-forward-loop motif features. The filename encodes the dump date; rename
  your download to this name or update `cfg.f_signor` in `run_all.m`.
- Reference: Licata et al., *NAR* 2020.

### `nci60_expr.txt` — CellMiner NCI-60 expression
- Source: https://discover.nci.nih.gov/cellminer (RNA-seq composite expression
  for the NCI-60 panel).
- Content: gene × cell-line expression matrix; cell-line naming matches ALMANAC.
- Reference: Reinhold et al., *Cancer Res.* 2012, 2019.

## Optional files

### `oneil_synergy.csv` — O'Neil screen (external validation, step 6)
- Source: DrugComb portal (study **ONEIL**), or the O'Neil subset of DrugCombDB.
- Expected columns (header, any order): `drug_row`, `drug_col`,
  `cell_line_name`, `synergy_zip`.
- If absent, `step6_external.m` prints a note and skips external validation.
- Reference: O'Neil et al., *Mol. Cancer Ther.* 2016.

### `drug_smiles.csv` and `drug_fingerprints.csv` — chemical fingerprints (steps 5, 7)
- `drug_smiles.csv`: two columns, `name,smiles`, with the pipeline's lower-case
  generic drug names (e.g. `dasatinib`, `5-fluorouracil`). Obtain SMILES from
  PubChem (https://pubchem.ncbi.nlm.nih.gov) or the DrugComb drug table.
- Generate `drug_fingerprints.csv` with:

  ```bash
  pip install rdkit
  python code/make_fingerprints.py      # reads data/drug_smiles.csv -> data/drug_fingerprints.csv
  ```

- If `drug_fingerprints.csv` is absent, the fingerprint baseline is skipped and
  the deep model falls back to target-profile features.
- References: Rogers & Hahn, *JCIM* 2010 (ECFP); Bento et al., *J. Cheminform.*
  2020 (RDKit).

## Note on column names

The loaders in `step1_features.m` and `step6_external.m` match a few columns by
header name. Database exports occasionally rename columns between releases; if a
step reports it cannot find a column, open the offending file, check the header,
and update the corresponding name in the step file (each is a short, clearly
commented block near the top of the function).
