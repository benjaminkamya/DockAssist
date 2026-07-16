"""
DockAssist — Ligand Pre-Screener for Molecular Docking

Accepts a ligand as:
1. Compound or drug name
2. SMILES
3. PubChem CID
4. PDB chemical component code

The input is converted to a molecular structure, analysed with RDKit,
and summarised using molecular descriptors and Lipinski-style checks.

Usage:
    python dock_assist.py
"""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any

import requests
from rdkit import Chem
from rdkit.Chem import Descriptors


PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
RCSB_BASE_URL = "https://data.rcsb.org/rest/v1/core/chemcomp"
REQUEST_TIMEOUT = 20


@dataclass(frozen=True)
class LigandInput:
    source: str
    original_value: str
    name: str
    smiles: str
    external_id: str | None = None


@dataclass(frozen=True)
class StructureMetadata:
    pdb_id: str
    title: str
    experimental_method: str
    resolution: float | None
    ligand_code: str


@dataclass(frozen=True)
class MolecularDescriptors:
    molecular_weight: float
    logp: float
    tpsa: float
    hydrogen_bond_donors: int
    hydrogen_bond_acceptors: int
    rotatable_bonds: int
    heavy_atoms: int
    ring_count: int


def request_json(url: str) -> dict[str, Any]:
    """Request JSON data and convert network problems into readable errors."""
    try:
        response = requests.get(url, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
    except requests.Timeout as error:
        raise ValueError("The external database request timed out.") from error
    except requests.ConnectionError as error:
        raise ValueError(
            "Could not connect to the external database. Check your internet connection."
        ) from error
    except requests.HTTPError as error:
        status = error.response.status_code if error.response is not None else "unknown"
        if status == 404:
            raise ValueError("No matching record was found.") from error
        raise ValueError(f"External database request failed with HTTP {status}.") from error

    try:
        return response.json()
    except ValueError as error:
        raise ValueError("The external database returned invalid JSON data.") from error


def validate_smiles(smiles: str) -> Chem.Mol:
    """Convert a SMILES string into an RDKit molecule."""
    cleaned = smiles.strip()

    if not cleaned:
        raise ValueError("SMILES input cannot be empty.")

    molecule = Chem.MolFromSmiles(cleaned)

    if molecule is None:
        raise ValueError("Invalid SMILES string.")

    return molecule


def canonicalise_smiles(smiles: str) -> str:
    """Return RDKit canonical SMILES for a valid structure."""
    molecule = validate_smiles(smiles)
    return Chem.MolToSmiles(molecule, canonical=True)


def pubchem_properties_from_url(url: str) -> dict[str, Any]:
    """Extract the first matching PubChem property record."""
    data = request_json(url)

    try:
        properties = data["PropertyTable"]["Properties"]
    except (KeyError, TypeError) as error:
        raise ValueError("PubChem returned an unexpected response.") from error

    if not properties:
        raise ValueError("No matching PubChem compound was found.")

    return properties[0]


def from_compound_name(value: str) -> LigandInput:
    """Resolve a compound name through PubChem."""
    encoded_name = requests.utils.quote(value.strip(), safe="")
    url = (
        f"{PUBCHEM_BASE_URL}/compound/name/{encoded_name}/property/"
        "Title,CanonicalSMILES,IsomericSMILES/JSON"
    )
    record = pubchem_properties_from_url(url)

    smiles = record.get("ConnectivitySMILES") or record.get(
        "CanonicalSMILES"
    ) or record.get("SMILES") or record.get("IsomericSMILES")

    if not smiles:
        raise ValueError("PubChem did not provide a usable SMILES string.")

    return LigandInput(
        source="PubChem compound name",
        original_value=value,
        name=record.get("Title") or value,
        smiles=canonicalise_smiles(smiles),
        external_id=str(record.get("CID")) if record.get("CID") else None,
    )


def from_pubchem_cid(value: str) -> LigandInput:
    """Resolve a PubChem CID."""
    try:
        cid = int(value.strip())
    except ValueError as error:
        raise ValueError("PubChem CID must be a whole number.") from error

    if cid <= 0:
        raise ValueError("PubChem CID must be greater than zero.")

    url = (
        f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/property/"
        "Title,CanonicalSMILES,IsomericSMILES/JSON"
    )
    record = pubchem_properties_from_url(url)

    smiles = record.get("ConnectivitySMILES") or record.get(
        "CanonicalSMILES"
    ) or record.get("SMILES") or record.get("IsomericSMILES")

    if not smiles:
        raise ValueError("PubChem did not provide a usable SMILES string.")

    return LigandInput(
        source="PubChem CID",
        original_value=value,
        name=record.get("Title") or f"PubChem CID {cid}",
        smiles=canonicalise_smiles(smiles),
        external_id=str(cid),
    )


def from_smiles(value: str) -> LigandInput:
    """Use a user-supplied SMILES string directly."""
    canonical = canonicalise_smiles(value)

    return LigandInput(
        source="SMILES",
        original_value=value,
        name="User-supplied molecule",
        smiles=canonical,
    )



def _is_likely_small_molecule_ligand(
    code: str,
    name: str,
    formula_weight: float | None,
) -> bool:
    """Exclude common solvents, ions and crystallisation additives."""
    excluded_codes = {
        "HOH", "DOD", "SO4", "PO4", "CL", "NA", "K", "CA", "MG", "MN",
        "ZN", "FE", "CU", "CO", "NI", "CD", "IOD", "BR", "GOL", "EDO",
        "PEG", "PG4", "ACT", "ACE", "FMT", "BME", "DMS", "MPD",
    }
    excluded_name_terms = {
        "water", "sulfate", "phosphate", "chloride", "sodium", "potassium",
        "calcium", "magnesium", "zinc", "glycerol", "ethylene glycol",
    }

    if code.upper() in excluded_codes:
        return False

    lowered_name = name.lower()
    if any(term in lowered_name for term in excluded_name_terms):
        return False

    if formula_weight is not None and formula_weight < 70:
        return False

    return True


def _extract_chemcomp_smiles(data: dict[str, Any]) -> str:
    """Extract the best available SMILES field from an RCSB chemcomp record."""
    descriptors = data.get("rcsb_chem_comp_descriptor", {})
    candidates = [
        descriptors.get("smiles_stereo"),
        descriptors.get("smiles"),
        descriptors.get("SMILES_stereo"),
        descriptors.get("SMILES"),
    ]

    for candidate in candidates:
        if candidate:
            return canonicalise_smiles(str(candidate))

    raise ValueError("RCSB did not provide a usable SMILES string for this ligand.")


def from_pdb_structure_id(
    value: str,
) -> tuple[LigandInput, StructureMetadata]:
    """Resolve a PDB structure and allow selection of a bound ligand."""
    pdb_id = value.strip().upper()

    if not re.fullmatch(r"[A-Z0-9]{4}", pdb_id):
        raise ValueError("PDB structure IDs contain exactly four letters or numbers.")

    entry = request_json(
        f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    )

    identifiers = entry.get("rcsb_entry_container_identifiers", {})
    entity_ids = identifiers.get("non_polymer_entity_ids") or []

    if not entity_ids:
        raise ValueError("No bound non-polymer ligands were found in this PDB entry.")

    ligand_options: list[dict[str, Any]] = []

    for entity_id in entity_ids:
        try:
            entity = request_json(
                "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/"
                f"{pdb_id}/{entity_id}"
            )
        except ValueError:
            continue

        id_data = entity.get(
            "rcsb_nonpolymer_entity_container_identifiers", {}
        )
        code = (
            id_data.get("nonpolymer_comp_id")
            or entity.get("pdbx_entity_nonpoly", {}).get("comp_id")
        )

        if not code:
            continue

        entity_data = entity.get("rcsb_nonpolymer_entity", {})
        name = (
            entity_data.get("pdbx_description")
            or entity.get("pdbx_entity_nonpoly", {}).get("name")
            or str(code)
        )

        formula_weight = entity_data.get("formula_weight")
        if formula_weight is not None:
            try:
                formula_weight = float(formula_weight) * 1000
            except (TypeError, ValueError):
                formula_weight = None

        if not _is_likely_small_molecule_ligand(
            str(code), str(name), formula_weight
        ):
            continue

        ligand_options.append(
            {
                "code": str(code).upper(),
                "name": str(name),
                "entity_id": str(entity_id),
            }
        )

    if not ligand_options:
        raise ValueError(
            "No likely small-molecule ligand could be identified in this entry."
        )

    if len(ligand_options) == 1:
        selected = ligand_options[0]
    else:
        print("\nBound ligands found:")
        for index, option in enumerate(ligand_options, start=1):
            print(f"{index}. {option['name']} ({option['code']})")

        selection = input("\nChoose ligand number: ").strip()

        try:
            selected = ligand_options[int(selection) - 1]
        except (ValueError, IndexError) as error:
            raise ValueError("Invalid ligand selection.") from error

    chemcomp = request_json(
        f"https://data.rcsb.org/rest/v1/core/chemcomp/{selected['code']}"
    )
    smiles = _extract_chemcomp_smiles(chemcomp)

    title = entry.get("struct", {}).get("title") or f"PDB structure {pdb_id}"

    methods = entry.get("exptl") or []
    experimental_method = ", ".join(
        str(item.get("method"))
        for item in methods
        if item.get("method")
    ) or "Not reported"

    resolution_values = (
        entry.get("rcsb_entry_info", {}).get("resolution_combined") or []
    )
    resolution = None
    if resolution_values:
        try:
            resolution = float(resolution_values[0])
        except (TypeError, ValueError):
            resolution = None

    ligand = LigandInput(
        source="RCSB PDB structure",
        original_value=value,
        name=selected["name"].title(),
        smiles=smiles,
        external_id=selected["code"],
    )

    structure = StructureMetadata(
        pdb_id=pdb_id,
        title=str(title),
        experimental_method=experimental_method,
        resolution=resolution,
        ligand_code=selected["code"],
    )

    return ligand, structure


def calculate_descriptors(molecule: Chem.Mol) -> MolecularDescriptors:
    """Calculate molecular descriptors relevant to ligand pre-screening."""
    return MolecularDescriptors(
        molecular_weight=Descriptors.MolWt(molecule),
        logp=Descriptors.MolLogP(molecule),
        tpsa=Descriptors.TPSA(molecule),
        hydrogen_bond_donors=Descriptors.NumHDonors(molecule),
        hydrogen_bond_acceptors=Descriptors.NumHAcceptors(molecule),
        rotatable_bonds=Descriptors.NumRotatableBonds(molecule),
        heavy_atoms=Descriptors.HeavyAtomCount(molecule),
        ring_count=Descriptors.RingCount(molecule),
    )


def check_lipinski(descriptors: MolecularDescriptors) -> dict[str, Any]:
    """Check the four conventional Lipinski thresholds."""
    checks = {
        "Molecular weight <= 500 Da": descriptors.molecular_weight <= 500,
        "LogP <= 5": descriptors.logp <= 5,
        "Hydrogen-bond donors <= 5": descriptors.hydrogen_bond_donors <= 5,
        "Hydrogen-bond acceptors <= 10": descriptors.hydrogen_bond_acceptors <= 10,
    }

    violations = sum(not passed for passed in checks.values())

    return {
        "checks": checks,
        "violations": violations,
        "passes": violations <= 1,
    }


def assess_ligand(descriptors: MolecularDescriptors) -> dict[str, Any]:
    """
    Perform a ligand-level, rule-based pre-screen.

    This does not predict receptor binding or docking performance.
    """
    checks = {
        "Molecular weight between 100 and 600 Da":
            100 <= descriptors.molecular_weight <= 600,
        "LogP between -1 and 6":
            -1 <= descriptors.logp <= 6,
        "TPSA <= 140 A^2":
            descriptors.tpsa <= 140,
        "Rotatable bonds <= 15":
            descriptors.rotatable_bonds <= 15,
        "At least 5 heavy atoms":
            descriptors.heavy_atoms >= 5,
    }

    failed = [description for description, passed in checks.items() if not passed]

    if not failed:
        recommendation = (
            "No obvious ligand-level physicochemical concerns were detected. "
            "Target-specific suitability still depends on the receptor and binding site."
        )
    elif len(failed) <= 2:
        recommendation = (
            "The ligand may still be dockable, but review the flagged properties before "
            "preparation. Target-specific suitability remains unknown."
        )
    else:
        recommendation = (
            "Several ligand-level properties may complicate preparation or docking. "
            "Manual review is recommended before continuing."
        )

    return {
        "checks": checks,
        "failed": failed,
        "recommendation": recommendation,
    }


def status_symbol(passed: bool) -> str:
    """Return a terminal-friendly status marker."""
    return "[PASS]" if passed else "[FAIL]"


def choose_input() -> tuple[LigandInput, StructureMetadata | None]:
    """Display the input menu and resolve the chosen ligand."""
    print("\nDockAssist")
    print("=" * 58)
    print("Choose an input type:")
    print("1. Compound or drug name")
    print("2. SMILES")
    print("3. PubChem CID")
    print("4. PDB structure ID")

    choice = input("\nSelection: ").strip()

    prompts = {
        "1": "Enter compound or drug name: ",
        "2": "Enter SMILES: ",
        "3": "Enter PubChem CID: ",
        "4": "Enter PDB structure ID: ",
    }

    if choice not in prompts:
        raise ValueError("Choose 1, 2, 3, or 4.")

    value = input(prompts[choice]).strip()

    resolvers = {
        "1": from_compound_name,
        "2": from_smiles,
        "3": from_pubchem_cid,
        "4": from_pdb_structure_id,
    }

    resolved = resolvers[choice](value)

    if choice == "4":
        return resolved

    return resolved, None


def display_results(
    ligand: LigandInput,
    structure: StructureMetadata | None,
    descriptors: MolecularDescriptors,
    lipinski_result: dict[str, Any],
    ligand_assessment: dict[str, Any],
) -> None:
    """Print the ligand pre-screen report."""
    print("\n" + "=" * 58)
    print("DockAssist — Ligand Pre-Screening Report")
    print("=" * 58)

    print("\nResolved Ligand")
    print("-" * 58)
    print(f"Name            : {ligand.name}")
    print(f"Input source    : {ligand.source}")
    print(f"Original input  : {ligand.original_value}")
    if ligand.external_id:
        print(f"External ID     : {ligand.external_id}")
    print(f"Canonical SMILES: {ligand.smiles}")

    if structure is not None:
        print("\nPDB Structure")
        print("-" * 58)
        print(f"PDB ID              : {structure.pdb_id}")
        print(f"Structure title     : {structure.title}")
        print(f"Experimental method : {structure.experimental_method}")
        if structure.resolution is not None:
            print(f"Resolution          : {structure.resolution:.2f} A")
        else:
            print("Resolution          : Not reported")
        print(f"Selected ligand     : {ligand.name} ({structure.ligand_code})")

    print("\nMolecular Descriptors")
    print("-" * 58)
    print(f"Molecular Weight : {descriptors.molecular_weight:.2f} Da")
    print(f"LogP             : {descriptors.logp:.2f}")
    print(f"TPSA             : {descriptors.tpsa:.2f} A^2")
    print(f"H-Bond Donors    : {descriptors.hydrogen_bond_donors}")
    print(f"H-Bond Acceptors : {descriptors.hydrogen_bond_acceptors}")
    print(f"Rotatable Bonds  : {descriptors.rotatable_bonds}")
    print(f"Heavy Atoms      : {descriptors.heavy_atoms}")
    print(f"Ring Count       : {descriptors.ring_count}")

    print("\nLipinski Assessment")
    print("-" * 58)
    for description, passed in lipinski_result["checks"].items():
        print(f"{status_symbol(passed):<8} {description}")

    print(f"\nViolations: {lipinski_result['violations']}")
    if lipinski_result["passes"]:
        print("Overall: Passes the conventional Lipinski assessment.")
    else:
        print("Overall: More than one Lipinski violation.")

    print("\nDocking-Oriented Ligand Assessment")
    print("-" * 58)
    for description, passed in ligand_assessment["checks"].items():
        print(f"{status_symbol(passed):<8} {description}")

    print("\nInterpretation")
    print("-" * 58)
    print(ligand_assessment["recommendation"])

    print("\nScientific Limitation")
    print("-" * 58)
    print(
        "DockAssist does not predict receptor binding, docking score, pose accuracy, "
        "selectivity, or biological activity."
    )


def main() -> int:
    """Run DockAssist."""
    try:
        ligand, structure = choose_input()
        molecule = validate_smiles(ligand.smiles)
        descriptors = calculate_descriptors(molecule)
        lipinski_result = check_lipinski(descriptors)
        ligand_assessment = assess_ligand(descriptors)

        display_results(
            ligand=ligand,
            structure=structure,
            descriptors=descriptors,
            lipinski_result=lipinski_result,
            ligand_assessment=ligand_assessment,
        )
    except ValueError as error:
        print(f"\nError: {error}")
        return 1
    except KeyboardInterrupt:
        print("\nOperation cancelled.")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
