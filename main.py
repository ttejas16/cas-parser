import pubchempy as pcp
import re
import json
from chempy import Substance
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from api import get_compound_info_from_cas


def get_equivalent_weight_rdkit(smiles_or_inchi, molecular_weight):
    if smiles_or_inchi.startswith("InChI"):
        mol = Chem.MolFromInchi(smiles_or_inchi)
    else:
        mol = Chem.MolFromSmiles(smiles_or_inchi)

    if mol:
        # Get formal charges to determine valency
        fragments = Chem.rdmolops.GetMolFrags(mol, asMols=True)
        max_charge = max(
            [abs(Chem.rdmolops.GetFormalCharge(frag)) for frag in fragments]
        )
        max_charge = max_charge if max_charge > 0 else 1
        return molecular_weight / max_charge


def get_equivalent_weight(formula, molecular_weight):
    substance = Substance.from_formula(formula)
    charge = substance.charge or 1
    return molecular_weight / abs(charge)


def get_cas_list():
    """reads cas numbers from a cas.txt"""

    with open("./cas.txt") as cas_nums:
        results = cas_nums.read().split("\n")
        return results


def get_simple_name(compound_info):
    if compound_info["synonyms"]:
        for syn in compound_info["synonyms"]:
            if re.match(r"^[a-zA-Z\s]+$", syn):
                return syn

    if compound_info["iupac_name"]:
        return compound_info["iupac_name"]

    # Last resort - CID
    return f"CID_{compound_info["cid"]}"


def is_salt_rdkit(smiles_or_inchi):
    """
    Determine if compound is a salt using RDKit structure analysis
    """
    # Convert to molecule object
    if smiles_or_inchi.startswith("InChI"):
        mol = Chem.MolFromInchi(smiles_or_inchi)
    else:
        mol = Chem.MolFromSmiles(smiles_or_inchi)

    if not mol:
        return False

    # Method 1: Fragment Analysis (most reliable)
    fragments = Chem.rdmolops.GetMolFrags(mol, asMols=True)

    if len(fragments) > 1:
        # Multiple fragments suggest ionic compound
        charged_fragments = 0
        cations = 0
        anions = 0

        for frag in fragments:
            charge = Chem.rdmolops.GetFormalCharge(frag)
            if charge != 0:
                charged_fragments += 1
                if charge > 0:
                    cations += 1
                else:
                    anions += 1

        # Salt if has both cations and anions
        if cations > 0 and anions > 0:
            return True

    # Method 2: Formal Charge Analysis
    total_formal_charge = Chem.rdmolops.GetFormalCharge(mol)
    atom_charges = [atom.GetFormalCharge() for atom in mol.GetAtoms()]

    # Check for charged atoms (zwitterions or ionic structures)
    positive_charges = sum(1 for charge in atom_charges if charge > 0)
    negative_charges = sum(1 for charge in atom_charges if charge < 0)

    if positive_charges > 0 and negative_charges > 0:
        return True

    return False


def is_salt_from_synonyms(compound_info):
    """
    Determine if compound is a salt by analyzing synonyms list
    Returns True if salt, False otherwise
    """
    if not compound_info["synonyms"]:
        return False

    # Convert all synonyms to lowercase for checking
    synonyms_lower = [syn.lower() for syn in compound_info["synonyms"]]
    synonyms_text = " ".join(synonyms_lower)

    # Strong salt indicators - these are definitive
    strong_salt_keywords = [
        "chloride",
        "bromide",
        "iodide",
        "fluoride",
        "sulfate",
        "sulphate",
        "nitrate",
        "phosphate",
        "carbonate",
        "bicarbonate",
        "hydroxide",
        "acetate",
        "citrate",
        "tartrate",
        "oxalate",
        "sodium",
        "potassium",
        "calcium",
        "magnesium",
        "lithium",
        "ammonium",
        "zinc",
        "copper",
        "hydrochloride",
        "hcl",
        "sulfonate",
        "benzoate",
        "succinate",
        "fumarate",
        "maleate",
        "lactate",
        "gluconate",
        "permanganate",
        "dichromate",
    ]

    # Check if any strong indicator is present
    for keyword in strong_salt_keywords:
        if keyword in synonyms_text:
            return True

    # Additional patterns that suggest salts
    salt_patterns = [
        "dihydrate",
        "trihydrate",
        "pentahydrate",
        "heptahydrate",  # Hydrated salts
        "monohydrate",
        "hexahydrate",
        "decahydrate",
        "anhydrous",  # Often used with salts
        "salt",
        "ionic",  # Direct mentions
    ]

    for pattern in salt_patterns:
        if pattern in synonyms_text:
            return True

    return False


result = {}
if __name__ == "__main__":
    cas_numbers = get_cas_list()
    not_found = 0
    with open("./skipped.txt", "w") as skipped:

        for cas in cas_numbers:
            compound = get_compound_info_from_cas(cas)

            if not compound:
                not_found += 1
                continue

            compound_name = get_simple_name(compound).lower()

            if not is_salt_from_synonyms(compound) or not is_salt_rdkit(
                compound["inchi"]
            ):
                skipped.write(f"skipped for cas: {cas} - {compound_name}\n")
                continue

            molecular_weight = float(compound["molecular_weight"])
            molecular_formula = compound["molecular_formula"]
            equivalent_weight = get_equivalent_weight_rdkit(
                compound["inchi"], molecular_weight
            )
            equivalent_weight = float(equivalent_weight)

            result[compound_name] = {
                "cas": cas,
                "molecular_weight": molecular_weight,
                "equivalent_weight": equivalent_weight,
                "molecular_formula": molecular_formula,
                "elements": compound["elements"],
                "iupac_name": compound["iupac_name"],
            }

    with open("result.json", "w") as f:
        json.dump(result, f)

    print(f"cannot search for {not_found} cas numbers")
