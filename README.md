# DockAssist

**DockAssist** is a Python tool for pre-screening small-molecule ligands before molecular docking. It accepts multiple forms of molecular input, calculates key physicochemical descriptors using RDKit, evaluates Lipinski's Rule of Five, and provides a ligand-level assessment to support early-stage docking workflows.

---

## Origin

DockAssist was inspired by my undergraduate dissertation 

The project began as an attempt to streamline ligand pre screening before molecular docking by automating descriptor calculation and providing rapid physicochemical assessment.

Rather than predicting binding affinity, DockAssist focuses on providing a rapid summary of molecular properties commonly examined during ligand selection. It is intended as a practical utility for computational drug discovery and structural biology workflows.
---

## Project Status

DockAssist is currently under active development.

The current version supports ligand retrieval through compound names, SMILES, PubChem CIDs and PDB Structure IDs. Molecular descriptors are calculated using RDKit before providing a ligand-level pre-screening assessment. Future versions will expand reporting, visualisation and structural analysis features.

---

## Features

- Accept ligand input through:
  - Compound or drug name
  - SMILES
  - PubChem CID
  - PDB Structure ID
- Retrieve ligand information from:
  - PubChem
  - RCSB Protein Data Bank
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
- Display PDB metadata including:
  - Structure ID
  - Structure title
  - Experimental method
  - Resolution
  - Bound ligand
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

Choose an input method:

```text
1. Compound Name
2. SMILES
3. PubChem CID
4. PDB Structure ID
```

---

## Example Output

```text
==========================================================
DockAssist
Ligand Pre-Screening Report
==========================================================

PDB ID              : 3PTB
Structure           : Trypsin
Experimental Method : X-Ray Diffraction
Resolution          : 1.80 Å

Ligand              : Benzamidine (BEN)

Molecular Weight    : 120.15 Da
LogP                : 0.97
TPSA                : 49.87 Å²
H-Bond Donors       : 2
H-Bond Acceptors    : 1
Rotatable Bonds     : 1

Lipinski Assessment
-------------------
PASS

Assessment
----------
No obvious ligand-level physicochemical concerns were detected.
```

---

## Roadmap

### Version 0.6.0 (Current)

- Compound name lookup
- SMILES input
- PubChem CID input
- PDB Structure ID support
- Protein metadata retrieval
- Molecular descriptor calculation
- Lipinski assessment

### Version 1.0.0

- Export reports
- 2D ligand visualisation
- Batch ligand screening

### Future Development

- Ligand similarity searching
- Protein-aware ligand comparison
- Machine learning-assisted ligand prioritisation

---

## Technologies

- Python
- RDKit
- PubChem
- RCSB Protein Data Bank

---

## Scientific Note

DockAssist is designed as a **ligand pre-screening tool**, not a docking engine. It evaluates molecular descriptors commonly considered before molecular docking studies but does **not** predict binding affinity, docking scores, binding poses, selectivity or biological activity.

---

## License

This project is released under the **MIT License**.
