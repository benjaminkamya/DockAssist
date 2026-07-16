"""
DockAssist. — Ligand Pre-Screener for Molecular Docking

Accepts a ligand as:
1. Compound or drug name
2. SMILES
3. PubChem CID
4. PDB chemical component code

The input is converted to a molecular structure, analysed with RDKit,
and summarised using molecular descriptors and Lipinski-style checks.

Usage:
    python ligand_scout.py
"""

from __future__ import annotations

from dataclasses import dataclass
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


def from_pdb_ligand_code(value: str) -> LigandInput:
    """Resolve an RCSB PDB chemical component code."""
    code = value.strip().upper()

    if not code:
        raise ValueError("PDB ligand code cannot be empty.")

    if len(code) > 5 or not code.isalnum():
        raise ValueError(
            "PDB ligand codes are short alphanumeric identifiers such as BEN or ATP."
        )

    data = request_json(f"{RCSB_BASE_URL}/{code}")

    descriptor_data = data.get("rcsb_chem_comp_descriptor", {})
    smiles = (
        descriptor_data.get("smiles_stereo")
        or descriptor_data.get("smiles")
    )

    if not smiles:
        raise ValueError("The PDB record did not provide a usable SMILES string.")

    chem_comp = data.get("chem_comp", {})
    name = chem_comp.get("name") or code

    return LigandInput(
        source="RCSB PDB chemical component",
        original_value=value,
        name=name.title(),
        smiles=canonicalise_smiles(smiles),
        external_id=code,
    )


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


def choose_input() -> LigandInput:
    """Display the input menu and resolve the chosen ligand."""
    print("\nLigandScout")
    print("=" * 58)
    print("Choose an input type:")
    print("1. Compound or drug name")
    print("2. SMILES")
    print("3. PubChem CID")
    print("4. PDB ligand code")

    choice = input("\nSelection: ").strip()

    prompts = {
        "1": "Enter compound or drug name: ",
        "2": "Enter SMILES: ",
        "3": "Enter PubChem CID: ",
        "4": "Enter PDB ligand code: ",
    }

    if choice not in prompts:
        raise ValueError("Choose 1, 2, 3, or 4.")

    value = input(prompts[choice]).strip()

    resolvers = {
        "1": from_compound_name,
        "2": from_smiles,
        "3": from_pubchem_cid,
        "4": from_pdb_ligand_code,
    }

    return resolvers[choice](value)


def display_results(
    ligand: LigandInput,
    descriptors: MolecularDescriptors,
    lipinski_result: dict[str, Any],
    ligand_assessment: dict[str, Any],
) -> None:
    """Print the ligand pre-screen report."""
    print("\n" + "=" * 58)
    print("LigandScout — Ligand Pre-Screening Report")
    print("=" * 58)

    print("\nResolved Ligand")
    print("-" * 58)
    print(f"Name            : {ligand.name}")
    print(f"Input source    : {ligand.source}")
    print(f"Original input  : {ligand.original_value}")
    if ligand.external_id:
        print(f"External ID     : {ligand.external_id}")
    print(f"Canonical SMILES: {ligand.smiles}")

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
        "LigandScout does not predict receptor binding, docking score, pose accuracy, "
        "selectivity, or biological activity."
    )


def main() -> int:
    """Run LigandScout."""
    try:
        ligand = choose_input()
        molecule = validate_smiles(ligand.smiles)
        descriptors = calculate_descriptors(molecule)
        lipinski_result = check_lipinski(descriptors)
        ligand_assessment = assess_ligand(descriptors)

        display_results(
            ligand=ligand,
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
