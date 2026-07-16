# DockAssist — Ligand Pre-Screening for Molecular Docking

DockAssist is a Python command-line tool that accepts a ligand by compound name, SMILES, PubChem CID, or PDB chemical component code. It resolves the molecular structure, calculates RDKit descriptors, applies Lipinski-style checks, and reports potential ligand-level concerns before molecular docking.

## Background

LigandScout was inspired by my undergraduate dissertation on molecular docking reproducibility. During that project, I found that there was no quick way to assess the physicochemical properties of ligands before beginning docking experiments.

This project was developed as a lightweight pre-screening tool to calculate molecular descriptors and summarise ligand properties prior to molecular docking. Future versions aim to incorporate structural information from the Protein Data Bank (PDB) to support protein–ligand workflows.

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
