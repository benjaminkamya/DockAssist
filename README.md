# DockAssist

**DockAssist** is a Python tool for pre-screening small-molecule ligands before molecular docking. It accepts multiple forms of molecular input, calculates key physicochemical descriptors using RDKit, evaluates Lipinski's Rule of Five, and provides a ligand-level assessment to support early-stage docking workflows.

---

## Background

DockAssist was inspired by my undergraduate dissertation investigating the reproducibility of molecular docking workflows. During that project, I often found myself wishing for a lightweight tool that could quickly characterise ligands before beginning docking experiments.

Rather than predicting binding affinity, DockAssist focuses on providing a rapid summary of molecular properties commonly examined during ligand selection. It is intended as a practical utility for computational drug discovery and structural biology workflows.

---

## Features

- Accept ligand input through:
  - Compound or drug name
  - SMILES
  - PubChem CID
  - *(Planned)* PDB Structure ID
- Calculate molecular descriptors using RDKit
  - Molecular Weight (MW)
  - LogP
  - Topological Polar Surface Area (TPSA)
  - Hydrogen Bond Donors (HBD)
  - Hydrogen Bond Acceptors (HBA)
  - Rotatable Bonds
  - Heavy Atom Count
  - Ring Count
- Evaluate Lipinski's Rule of Five
- Provide a rule-based ligand pre-screening assessment
- Simple command-line interface

---

## Installation

Install the required packages:

```bash
pip install -r requirements.txt
```

---

## Usage

Run the program:

```bash
python dock_assist.py
```

Choose an input method and enter the requested information.

Example:

```
1. Compound Name
2. SMILES
3. PubChem CID
4. PDB Structure ID (Coming Soon)
```

---

## Example Output

```
==========================================================
DockAssist
Ligand Pre-Screening Report
==========================================================

Name            : Benzamidine

Molecular Weight : 120.15 Da
LogP             : 0.97
TPSA             : 49.87 Å²
H-Bond Donors    : 2
H-Bond Acceptors : 1
Rotatable Bonds  : 1

Lipinski Assessment
-------------------
PASS

Assessment
----------
No obvious ligand-level physicochemical concerns were detected.

Scientific Limitation
---------------------
DockAssist does not predict receptor binding, docking scores,
binding affinity, selectivity or biological activity.
```

---

## Roadmap

### Version 0.5.0 (Current)

- Compound name lookup
- SMILES input
- PubChem CID input
- Molecular descriptor calculation
- Lipinski assessment

### Version 1.0.0

- PDB Structure ID support
- Protein metadata retrieval
- Automatic identification of bound ligands
- Ligand descriptor analysis from PDB structures

### Future Development

- 2D ligand visualisation
- CSV export
- Batch ligand screening
- Ligand similarity searching
- Protein-aware ligand comparison
- Machine learning-assisted ligand prioritisation

---

## Technologies

- Python
- RDKit
- PubChem
- RCSB Protein Data Bank (planned)

---

## Scientific Note

DockAssist is designed as a **ligand pre-screening tool**, not a docking engine. It evaluates molecular descriptors that are commonly considered before molecular docking studies but does **not** predict binding affinity or docking outcomes.

---

## License

MIT License
