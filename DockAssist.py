"""
DockAssist — Ligand Pre-Screener for Molecular Docking

Screens ligands with RDKit descriptors, Lipinski's rules and docking-oriented
checks. Input can be a compound name, SMILES, PubChem CID or PDB structure ID.

Usage:
    python DockAssist.py                                 # menu (single or batch)
    python DockAssist.py --batch ligands.csv [--out dir] # batch straight away

In batch mode the chemical-space plot is saved to the output folder and opened.

Batch CSV columns: "input" (required: name, SMILES or CID), "type" (optional:
name, smiles or cid; guessed if blank) and "name" (optional display label).
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
import os
from pathlib import Path
import re
import subprocess
import sys
import time

import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")  # bad SMILES are reported by DockAssist itself

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
RCSB = "https://data.rcsb.org/rest/v1/core"
TIMEOUT = 20
PUBCHEM_PAUSE = 0.25  # PubChem allows at most 5 requests per second

# Common solvents, ions and crystallisation additives to ignore in PDB entries
NON_LIGAND_CODES = {
    "HOH", "DOD", "SO4", "PO4", "CL", "NA", "K", "CA", "MG", "MN", "ZN", "FE",
    "CU", "CO", "NI", "CD", "IOD", "BR", "GOL", "EDO", "PEG", "PG4", "ACT",
    "ACE", "FMT", "BME", "DMS", "MPD",
}
NON_LIGAND_TERMS = (
    "water", "sulfate", "phosphate", "chloride", "sodium", "potassium",
    "calcium", "magnesium", "zinc", "glycerol", "ethylene glycol",
)

DESCRIPTORS = {
    "molecular_weight": Descriptors.MolWt,
    "logp": Descriptors.MolLogP,
    "tpsa": Descriptors.TPSA,
    "hbd": Descriptors.NumHDonors,
    "hba": Descriptors.NumHAcceptors,
    "rotatable_bonds": Descriptors.NumRotatableBonds,
    "heavy_atoms": Descriptors.HeavyAtomCount,
    "ring_count": Descriptors.RingCount,
}

LIPINSKI_RULES = {
    "Molecular weight <= 500 Da": lambda d: d["molecular_weight"] <= 500,
    "LogP <= 5": lambda d: d["logp"] <= 5,
    "Hydrogen-bond donors <= 5": lambda d: d["hbd"] <= 5,
    "Hydrogen-bond acceptors <= 10": lambda d: d["hba"] <= 10,
}

DOCKING_RULES = {
    "Molecular weight between 100 and 600 Da": lambda d: 100 <= d["molecular_weight"] <= 600,
    "LogP between -1 and 6": lambda d: -1 <= d["logp"] <= 6,
    "TPSA <= 140 A^2": lambda d: d["tpsa"] <= 140,
    "Rotatable bonds <= 15": lambda d: d["rotatable_bonds"] <= 15,
    "At least 5 heavy atoms": lambda d: d["heavy_atoms"] >= 5,
}


@dataclass(frozen=True)
class Ligand:
    name: str
    source: str
    original_value: str
    smiles: str
    external_id: str | None = None


@dataclass(frozen=True)
class Structure:
    pdb_id: str
    title: str
    method: str
    resolution: float | None
    ligand_code: str


# ---------------------------------------------------------------------------
# Resolving inputs
# ---------------------------------------------------------------------------

def get_json(url: str) -> dict:
    """GET a URL and turn network problems into readable errors."""
    try:
        response = requests.get(url, timeout=TIMEOUT)
        response.raise_for_status()
        return response.json()
    except requests.Timeout as error:
        raise ValueError("The external database request timed out.") from error
    except requests.ConnectionError as error:
        raise ValueError("Could not connect to the external database. "
                         "Check your internet connection.") from error
    except requests.HTTPError as error:
        status = getattr(error.response, "status_code", "unknown")
        if status == 404:
            raise ValueError("No matching record was found.") from error
        raise ValueError(f"External database request failed with HTTP {status}.") from error
    except ValueError as error:
        raise ValueError("The external database returned invalid JSON data.") from error


def to_mol(smiles: str) -> Chem.Mol:
    """Parse SMILES into an RDKit molecule, or raise a readable error."""
    molecule = Chem.MolFromSmiles(smiles.strip()) if smiles.strip() else None
    if molecule is None:
        raise ValueError("Invalid or empty SMILES string.")
    return molecule


def canonical(smiles: str) -> str:
    return Chem.MolToSmiles(to_mol(smiles))


def from_smiles(value: str) -> Ligand:
    return Ligand("User-supplied molecule", "SMILES", value, canonical(value))


def from_pubchem(value: str, by: str) -> Ligand:
    """Look up a compound in PubChem by 'name' or 'cid'."""
    if by == "cid":
        if not value.strip().isdigit() or int(value) <= 0:
            raise ValueError("PubChem CID must be a positive whole number.")
        key = value.strip()
    else:
        key = requests.utils.quote(value.strip(), safe="")

    data = get_json(f"{PUBCHEM}/{by}/{key}/property/"
                    "Title,CanonicalSMILES,IsomericSMILES/JSON")
    try:
        record = data["PropertyTable"]["Properties"][0]
    except (KeyError, IndexError, TypeError) as error:
        raise ValueError("No matching PubChem compound was found.") from error

    # PubChem has renamed its SMILES fields over time, so try each in turn
    smiles = next((record[k] for k in ("ConnectivitySMILES", "CanonicalSMILES",
                                       "SMILES", "IsomericSMILES") if record.get(k)), None)
    if not smiles:
        raise ValueError("PubChem did not provide a usable SMILES string.")

    cid = str(record.get("CID") or value.strip())
    source = "PubChem CID" if by == "cid" else "PubChem compound name"
    return Ligand(record.get("Title") or value, source, value, canonical(smiles), cid)


def from_pdb(value: str) -> tuple[Ligand, Structure]:
    """Resolve a PDB entry and let the user pick a bound ligand."""
    pdb_id = value.strip().upper()
    if not re.fullmatch(r"[A-Z0-9]{4}", pdb_id):
        raise ValueError("PDB structure IDs contain exactly four letters or numbers.")

    entry = get_json(f"{RCSB}/entry/{pdb_id}")
    entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get(
        "non_polymer_entity_ids") or []

    options = []
    for entity_id in entity_ids:
        try:
            entity = get_json(f"{RCSB}/nonpolymer_entity/{pdb_id}/{entity_id}")
        except ValueError:
            continue
        nonpoly = entity.get("pdbx_entity_nonpoly", {})
        info = entity.get("rcsb_nonpolymer_entity", {})
        code = (entity.get("rcsb_nonpolymer_entity_container_identifiers", {})
                .get("nonpolymer_comp_id") or nonpoly.get("comp_id"))
        if not code:
            continue
        name = str(info.get("pdbx_description") or nonpoly.get("name") or code)
        try:
            weight = float(info["formula_weight"]) * 1000  # stored in kDa
        except (KeyError, TypeError, ValueError):
            weight = None

        is_additive = (code.upper() in NON_LIGAND_CODES
                       or any(term in name.lower() for term in NON_LIGAND_TERMS)
                       or (weight is not None and weight < 70))
        if not is_additive:
            options.append((code.upper(), name))

    if not options:
        raise ValueError("No likely small-molecule ligand was found in this entry.")

    if len(options) == 1:
        code, name = options[0]
    else:
        print("\nBound ligands found:")
        for i, (code, name) in enumerate(options, start=1):
            print(f"{i}. {name} ({code})")
        try:
            code, name = options[int(input("\nChoose ligand number: ")) - 1]
        except (ValueError, IndexError) as error:
            raise ValueError("Invalid ligand selection.") from error

    descriptors = get_json(f"{RCSB}/chemcomp/{code}").get("rcsb_chem_comp_descriptor", {})
    smiles = next((descriptors[k] for k in ("smiles_stereo", "smiles", "SMILES_stereo",
                                            "SMILES") if descriptors.get(k)), None)
    if not smiles:
        raise ValueError("RCSB did not provide a usable SMILES string for this ligand.")

    resolution = (entry.get("rcsb_entry_info", {}).get("resolution_combined") or [None])[0]
    methods = ", ".join(m["method"] for m in entry.get("exptl") or [] if m.get("method"))

    ligand = Ligand(name.title(), "RCSB PDB structure", value, canonical(smiles), code)
    structure = Structure(
        pdb_id=pdb_id,
        title=entry.get("struct", {}).get("title") or f"PDB structure {pdb_id}",
        method=methods or "Not reported",
        resolution=float(resolution) if resolution is not None else None,
        ligand_code=code,
    )
    return ligand, structure


RESOLVERS = {
    "name": lambda v: from_pubchem(v, "name"),
    "smiles": from_smiles,
    "cid": lambda v: from_pubchem(v, "cid"),
}


def guess_type(value: str) -> str:
    """Whole numbers are CIDs, parseable strings are SMILES, the rest are names."""
    if value.strip().isdigit():
        return "cid"
    return "smiles" if Chem.MolFromSmiles(value.strip()) else "name"


# ---------------------------------------------------------------------------
# Screening
# ---------------------------------------------------------------------------

def screen(ligand: Ligand) -> dict:
    """Calculate descriptors and rule checks for one ligand as a flat row."""
    molecule = to_mol(ligand.smiles)
    d = {key: round(func(molecule), 2) for key, func in DESCRIPTORS.items()}
    violations = sum(not rule(d) for rule in LIPINSKI_RULES.values())
    flags = [label for label, rule in DOCKING_RULES.items() if not rule(d)]
    return {
        "name": ligand.name, "source": ligand.source,
        "external_id": ligand.external_id, "smiles": ligand.smiles, **d,
        "lipinski_violations": violations,
        "passes_lipinski": violations <= 1,
        "docking_flags": len(flags),
        "flagged_properties": "; ".join(flags),
    }


def print_report(ligand: Ligand, row: dict, structure: Structure | None) -> None:
    def section(title: str) -> None:
        print(f"\n{title}\n" + "-" * 58)

    def checks(rules: dict) -> None:
        for label, rule in rules.items():
            print(f"{'[PASS]' if rule(row) else '[FAIL]':<8} {label}")

    print("\n" + "=" * 58 + "\nDockAssist — Ligand Pre-Screening Report\n" + "=" * 58)

    section("Resolved Ligand")
    print(f"Name            : {ligand.name}\nInput source    : {ligand.source}\n"
          f"Original input  : {ligand.original_value}")
    if ligand.external_id:
        print(f"External ID     : {ligand.external_id}")
    print(f"Canonical SMILES: {ligand.smiles}")

    if structure:
        section("PDB Structure")
        res = f"{structure.resolution:.2f} A" if structure.resolution else "Not reported"
        print(f"PDB ID              : {structure.pdb_id}\n"
              f"Structure title     : {structure.title}\n"
              f"Experimental method : {structure.method}\n"
              f"Resolution          : {res}\n"
              f"Selected ligand     : {ligand.name} ({structure.ligand_code})")

    section("Molecular Descriptors")
    print(f"Molecular Weight : {row['molecular_weight']:.2f} Da\n"
          f"LogP             : {row['logp']:.2f}\n"
          f"TPSA             : {row['tpsa']:.2f} A^2\n"
          f"H-Bond Donors    : {row['hbd']}\nH-Bond Acceptors : {row['hba']}\n"
          f"Rotatable Bonds  : {row['rotatable_bonds']}\n"
          f"Heavy Atoms      : {row['heavy_atoms']}\nRing Count       : {row['ring_count']}")

    section("Lipinski Assessment")
    checks(LIPINSKI_RULES)
    print(f"\nViolations: {row['lipinski_violations']}\nOverall: " +
          ("Passes the conventional Lipinski assessment." if row["passes_lipinski"]
           else "More than one Lipinski violation."))

    section("Docking-Oriented Ligand Assessment")
    checks(DOCKING_RULES)

    section("Interpretation")
    n = row["docking_flags"]
    if n == 0:
        print("No obvious ligand-level physicochemical concerns were detected. "
              "Target-specific suitability still depends on the receptor and binding site.")
    elif n <= 2:
        print("The ligand may still be dockable, but review the flagged properties "
              "before preparation. Target-specific suitability remains unknown.")
    else:
        print("Several ligand-level properties may complicate preparation or docking. "
              "Manual review is recommended before continuing.")

    section("Scientific Limitation")
    print("DockAssist does not predict receptor binding, docking score, pose accuracy, "
          "selectivity, or biological activity.")


def run_interactive() -> int:
    print("\nDockAssist\n" + "=" * 58 + "\nChoose an input type:\n"
          "1. Compound or drug name\n2. SMILES\n3. PubChem CID\n4. PDB structure ID\n"
          "5. Batch screen a CSV file (with chemical-space plot)")
    menu = {"1": "name", "2": "smiles", "3": "cid", "4": "pdb", "5": "batch"}
    try:
        choice = menu.get(input("\nSelection: ").strip())
        if choice is None:
            raise ValueError("Choose 1, 2, 3, 4 or 5.")
        if choice == "batch":
            default = Path(__file__).with_name("example_ligands.csv")
            entered = input(f"CSV file [Enter for {default.name}]: ").strip().strip('"')
            csv_path = Path(entered) if entered else default
            return run_batch(csv_path, csv_path.parent / "results")
        value = input(f"Enter {'PDB structure ID' if choice == 'pdb' else choice}: ").strip()
        ligand, structure = from_pdb(value) if choice == "pdb" else (RESOLVERS[choice](value), None)
        print_report(ligand, screen(ligand), structure)
    except ValueError as error:
        print(f"\nError: {error}")
        return 1
    except KeyboardInterrupt:
        print("\nOperation cancelled.")
        return 1
    return 0


# ---------------------------------------------------------------------------
# Batch mode
# ---------------------------------------------------------------------------

def run_batch(csv_path: Path, out_dir: Path) -> int:
    """Screen every ligand in a CSV, clean the results, save table and plot."""
    import pandas as pd

    try:
        table = pd.read_csv(csv_path, dtype=str).fillna("")
    except FileNotFoundError:
        print(f"Error: could not find {csv_path}")
        return 1
    if "input" not in table.columns:
        print('Error: the CSV needs a column called "input".')
        return 1
    for column in ("type", "name"):
        table[column] = table.get(column, "")

    # Cleaning 1: drop blank rows before any lookups
    table["input"] = table["input"].str.strip()
    table = table[table["input"] != ""]

    rows, errors = [], []
    print(f"\nScreening {len(table)} ligands from {csv_path.name}...")
    for record in table.itertuples():
        kind = record.type.strip().lower() or guess_type(record.input)
        try:
            if kind not in RESOLVERS:
                raise ValueError(f"Unknown type '{kind}'.")
            ligand = RESOLVERS[kind](record.input)
            if record.name.strip():
                ligand = replace(ligand, name=record.name.strip())
            rows.append({"input": record.input, **screen(ligand)})
            print(f"  [OK]   {ligand.name}")
        except ValueError as error:
            # Cleaning 2: a failed lookup or bad SMILES is logged, not fatal
            errors.append({"input": record.input, "error": str(error)})
            print(f"  [SKIP] {record.input}: {error}")
        if kind in ("name", "cid"):
            time.sleep(PUBCHEM_PAUSE)

    if not rows:
        print("\nNo ligands could be screened.")
        return 1

    # Cleaning 3: the same molecule entered twice has the same canonical SMILES
    results = pd.DataFrame(rows)
    duplicates = results.duplicated(subset="smiles").sum()
    results = results.drop_duplicates(subset="smiles")

    out_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(out_dir / "dockassist_results.csv", index=False)
    if errors:
        pd.DataFrame(errors).to_csv(out_dir / "dockassist_errors.csv", index=False)
    plot_path = out_dir / "dockassist_chemical_space.png"
    plot_chemical_space(results, plot_path)

    n = len(results)
    medians = results[["molecular_weight", "logp", "tpsa", "rotatable_bonds"]].median()
    print("\n" + "=" * 58 + "\nDockAssist — Batch Summary\n" + "=" * 58)
    print(f"Ligands in file           : {len(table)}\n"
          f"Screened                  : {n}\n"
          f"Skipped (errors)          : {len(errors)}\n"
          f"Duplicates removed        : {duplicates}\n"
          f"Pass Lipinski (<=1 viol.) : {results['passes_lipinski'].sum()} of {n}\n"
          f"No docking-oriented flags : {(results['docking_flags'] == 0).sum()} of {n}")
    print("\nMedian descriptors of screened set\n" + "-" * 58)
    print(medians.round(2).to_string())
    print(f"\nOutputs saved to {out_dir.resolve()}")
    open_file(plot_path)
    return 0


def open_file(path: Path) -> None:
    """Open a file in the system's default viewer (best effort)."""
    try:
        if sys.platform.startswith("win"):
            os.startfile(path)  # type: ignore[attr-defined]
        elif sys.platform == "darwin":
            subprocess.run(["open", str(path)], check=False)
        else:
            subprocess.run(["xdg-open", str(path)], check=False,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except OSError:
        print(f"Open the plot here: {path.resolve()}")


def plot_chemical_space(results, plot_path: Path) -> None:
    """Scatter MW against logP with the Lipinski region marked."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, muted = "#fcfcfb", "#1f1e1d", "#6b6a66"
    blue, orange = "#2a78d6", "#eb6834"

    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=200, facecolor=surface)
    ax.set_facecolor(surface)
    y_min = min(-2.0, results["logp"].min() - 0.8)
    y_max = max(7.0, results["logp"].max() + 0.8)

    ax.fill_between([0, 500], y_min, 5, color=blue, alpha=0.07, linewidth=0)
    ax.axvline(500, color=muted, linewidth=1, linestyle="--")
    ax.axhline(5, color=muted, linewidth=1, linestyle="--")
    ax.text(15, 4.75, "MW ≤ 500 Da and logP ≤ 5", color=muted, fontsize=9, va="top")

    passed = results[results["passes_lipinski"]]
    failed = results[~results["passes_lipinski"]]
    for subset, colour, marker, label in (
        (passed, blue, "o", f"Passes Lipinski, 0–1 violations ({len(passed)})"),
        (failed, orange, "D", f"2+ Lipinski violations ({len(failed)})"),
    ):
        ax.scatter(subset["molecular_weight"], subset["logp"], s=75, marker=marker,
                   color=colour, edgecolor=surface, linewidth=1.5, label=label, zorder=3)

    for row in failed.itertuples():  # label only the failures
        ax.annotate(row.name, (row.molecular_weight, row.logp), xytext=(7, 5),
                    textcoords="offset points", fontsize=9, color=ink)

    ax.set(xlim=(0, max(650.0, results["molecular_weight"].max() * 1.08)),
           ylim=(y_min, y_max))
    ax.set_xlabel("Molecular weight (Da)", color=ink)
    ax.set_ylabel("logP (Crippen)", color=ink)
    ax.set_title("DockAssist: screened ligands in chemical space",
                 color=ink, fontsize=13, loc="left", pad=12)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(muted)
    ax.tick_params(colors=muted)
    ax.grid(color=muted, alpha=0.15, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", frameon=False, fontsize=9, labelcolor=ink)

    fig.tight_layout()
    fig.savefig(plot_path, facecolor=surface)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="DockAssist ligand pre-screener")
    parser.add_argument("--batch", type=Path, help="CSV of ligands ('input' column)")
    parser.add_argument("--out", type=Path, default=Path("results"),
                        help="folder for batch outputs (default: results)")
    args = parser.parse_args()
    return run_batch(args.batch, args.out) if args.batch else run_interactive()


if __name__ == "__main__":
    raise SystemExit(main())
