#!/usr/bin/env python3
"""
make_fingerprints.py  -  turn a drug->SMILES table into Morgan fingerprints.

Use this only if you want the chemical-fingerprint bar in Figure 5. It needs
RDKit (pip install rdkit).

Input : data/drug_smiles.csv   with a header and two columns: name,smiles
Output: data/drug_fingerprints.csv   name followed by the fingerprint bits

The drug names must match the pipeline's names (lower-case generics, e.g.
"dasatinib", "5-fluorouracil"); step5_baseline.m matches on them.

Get SMILES from any source your browser can reach, e.g. PubChem
(pubchem.ncbi.nlm.nih.gov, "Download > CSV") or the DrugComb drug table.
"""
import csv, sys, os
from rdkit import Chem
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), 'data')
INP  = os.path.join(DATA, 'drug_smiles.csv')
OUT  = os.path.join(DATA, 'drug_fingerprints.csv')
NBITS = 1024
RADIUS = 2

def main():
    if not os.path.exists(INP):
        sys.exit(f'missing {INP} (two columns: name,smiles)')
    rows, bad = [], []
    with open(INP, newline='') as fh:
        r = csv.reader(fh); next(r, None)
        for line in r:
            if len(line) < 2: continue
            name, smi = line[0].strip().lower(), line[1].strip()
            m = Chem.MolFromSmiles(smi)
            if m is None:
                bad.append(name); continue
            fp = AllChem.GetMorganFingerprintAsBitVect(m, RADIUS, nBits=NBITS)
            rows.append([name] + list(fp))
    with open(OUT, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['name'] + [f'b{i}' for i in range(NBITS)])
        w.writerows(rows)
    print(f'wrote {len(rows)} fingerprints -> {OUT}')
    if bad:
        print(f'could not parse {len(bad)}: {", ".join(bad)}')

if __name__ == '__main__':
    main()
