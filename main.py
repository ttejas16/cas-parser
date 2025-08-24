import pubchempy as pcp
import re
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


if __name__ == "__main__":
    cas_numbers = get_cas_list()

    with open("./result.txt", "w") as result:
        result.write("cas,name,molecular_weight,equivalent_weight\n")

        for cas in cas_numbers:
            compound = get_compound_info_from_cas(cas)

            if not compound:
                print("Could not find for cas:", cas)
                continue

            compound_name = get_simple_name(compound).lower()
            molecular_weight = float(compound["molecular_weight"])
            molecular_formula = compound["molecular_formula"]
            equivalent_weight = get_equivalent_weight_rdkit(
                compound["inchi"], molecular_weight
            )

            writable_string = "\t\t".join(
                map(str, [cas, compound_name, molecular_weight, equivalent_weight])
            )
            result.write(writable_string + "\n")
