import pubchempy as pcp
import re
from chempy import Substance
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


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


def get_simple_name(compound):
    if compound.synonyms:
        for syn in compound.synonyms:
            if re.match(r'^[a-zA-Z\s]+$', syn):
                return syn
    
    if compound.iupac_name:
        return compound.iupac_name

    # Last resort - CID
    return f"CID_{compound.cid}"


if __name__ == "__main__":
    cas_numbers = get_cas_list()

    with open("./result.txt", "w") as result:
        result.write("cas,name,molecular_weight,equivalent_weight\n")

        for cas in cas_numbers:
            compounds = pcp.get_compounds(cas, "name")
            if not compounds:
                print("Could not find for cas:", cas)
                continue

            compound = compounds[0]

            compound_name = get_simple_name(compound).lower()
            molecular_weight = float(compound.molecular_weight)
            molecular_formula = compound.molecular_formula
            equivalent_weight = get_equivalent_weight_rdkit(
                compound.inchi, molecular_weight
            )

            writable_string = "\t\t".join(
                map(str, [cas, compound_name, molecular_weight, equivalent_weight])
            )
            result.write(writable_string + "\n")
