from mp_api.client import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import json

molecule_species_mpid = 'mp-22893'

with MPRester('API-KEY') as mpr:
    primitive_structure = mpr.get_structure_by_material_id(molecule_species_mpid)

#conv_struct = SpacegroupAnalyzer(primitive_structure).get_conventional_standard_structure()

with open('mp22893_primitive.json','w') as f:
    json.dump(primitive_structure.as_dict(), f)

