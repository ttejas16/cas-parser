import pubchempy as pcp
import os
import json

cache_file_name = ".cache.json"
cache = None

# writing this because cache should always exists with an empty object on startup
if not os.path.exists(cache_file_name):
    with open(cache_file_name, "w") as f:
        f.write("{}")
else:
    with open(cache_file_name, "r+") as f:
        if len(f.readline()) == 0:
            f.write("{}")


with open(cache_file_name, "r") as cf:
    cache = json.load(cf)



def save_to_cache(cas, info):
    """save the cas and its info to local cache and update global cache object"""
    global cache
    cache[cas] = info

    with open(cache_file_name, "w") as f:
        json.dump(cache, f)


def get_compound_info_from_cas(cas):
    """
    Caches and retrieve compound information for CAS numbers.

    Args:
        cas (str): CAS number in format 'XXXXX-XX-X'

    Returns:
        dict: Compound information with the following structure:
            {
                'molecular_weight': str,
                'synonyms': List[str],
                'inchi': str,
                'molecular_formula': str,
                'cid': str,
                'iupac_name': str,
                'elements': List[str]
            }
        Returns None if compound not found.
    """

    entry = cache.get(cas, None)
    if entry:
        # print("Cache hit for cas: ", cas)
        return entry

    compounds = pcp.get_compounds(cas, "name")
    if not compounds:
        return None

    compound = compounds[0]

    obj = {}
    obj["cid"] = compound.cid
    obj["synonyms"] = compound.synonyms
    obj["molecular_weight"] = compound.molecular_weight
    obj["molecular_formula"] = compound.molecular_formula
    obj["inchi"] = compound.inchi
    obj["iupac_name"] = compound.iupac_name
    obj["elements"] = compound.elements

    save_to_cache(cas, obj)
    return obj
