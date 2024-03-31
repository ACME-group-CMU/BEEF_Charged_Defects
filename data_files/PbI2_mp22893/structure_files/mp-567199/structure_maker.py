from mp_api.client import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import json

molecule_species_mpid = 'mp-567199'
primitive = True

with MPRester('API_KEY') as mpr:
    primitive_structure = mpr.get_structure_by_material_id(molecule_species_mpid)

#conv_struct = SpacegroupAnalyzer(primitive_structure).get_conventional_standard_structure()

if(primitive):
	with open(f"{molecule_species_mpid}_primitive",'w') as f:
    		json.dump(primitive_structure.as_dict(), f)
else:
	conv_struct = SpacegroupAnalyzer(primitive_structure).get_conventional_standard_structure()
	with open(f"{molecule_species_mpid}_conventional",'w') as f:
    		json.dump(conv_struct.as_dict(), f)

