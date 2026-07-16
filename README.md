# DockAssist — Ligand Pre-Screening for Molecular Docking

DockAssist is a Python command-line tool that accepts a ligand by compound name, SMILES, PubChem CID, or PDB chemical component code. It resolves the molecular structure, calculates RDKit descriptors, applies Lipinski-style checks, and reports potential ligand-level concerns before molecular docking.

## Input methods

- Compound or drug name through PubChem
- SMILES
- PubChem CID
- PDB ligand code through the RCSB Protein Data Bank

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python ligand_scout.py
```

Example PDB ligand input:

```text
Selection: 4
Enter PDB ligand code: BEN
```

Example compound-name input:

```text
Selection: 1
Enter compound or drug name: benzamidine
```

## Scientific limitation

DockAssist performs ligand-level physicochemical pre-screening. It does not predict receptor binding, docking score, pose accuracy, selectivity, or biological activity.
