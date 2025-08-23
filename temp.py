import pubchempy as pcp

compound = pcp.get_compounds("6381-92-6", 'name')[0]

print(compound.iupac_name)
print(compound.synonyms)