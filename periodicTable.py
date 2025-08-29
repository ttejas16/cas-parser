import pubchempy as pcp
import re
from pubchempy import ELEMENTS
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from api import get_element_info_from_name
import time

allElements = list(ELEMENTS.values())
data = {}
for e in allElements:
    print("Fetching data for element:", e)
    
    data[e] = get_element_info_from_name(e)


import json

# Write the data dictionary to a JSON file
with open("element_data.json", "w", encoding="utf-8") as f:
    json.dump(data, f, ensure_ascii=False, indent=4)