#Pymatgen
from pymatgen.analysis.defects.generators import VacancyGenerator
from pymatgen.ext.matproj import MPRester

#Other
import importlib
import os
import subprocess
import shutil

#Defect Tracking
from qe_defect_tracker import Tracker
importlib.reload(Tracker)

os.environ['ESPRESSO_PW_EXE'] = "pw.x"
os.environ['ESPRESSO_PP_EXE'] = "pp.x"
os.environ['COFFEE_DIR'] = "/trace/home/atimmins/packages/coffee/CoFFEE-main/"

#Define relevant details
molecule_species_mpid = 'mp-48'
force_recalc = False

#K-points (can check multiple)
kpts = [10,10,10] #update

# Define pseudo potentials
pseudopotentials = {'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF'}

#Generate the ideal dimensions of the supercell (pymatgen uses this ASE function, could incorporate more cleanly in the future
#conf = bulk(molecule_species)
#supercell_dim = find_optimal_cell_shape(conf.cell, molecule_number, 'sc')
supercell_dim = [1,1,1]

# Define the input data specific to Quantum Espresso
input_data = {
    'control' : {
        'calculation': 'ensemble',
        'outdir':'./',
        'pseudo_dir': '/trace/group/acmegroup/atimmins/packages/espresso/pseudo/psl_1.0.0/pslibrary.1.0.0/pbe/PSEUDOPOTENTIALS/',
        'prefix':'GaN',
        'verbosity':'high',
        'restart_mode': 'from_scratch',
	'tstress': True,
	'tprnfor': True,
    },
    'system' : {
        'nat': 1,
        'ntyp': 1,
        'ecutwfc':70,
	'ecutrho':350,
        'tot_charge':0,
        'ibrav': 4,
	'degauss':0.01,
	'smearing':'gauss',
	'occupations':'smearing',
	'input_dft':'BEEF-vdW',
    },
    'electrons': {
	'electron_maxstep':100,
	'conv_thr': 1e-8,
	'mixing_mode': 'plain',
	'mixing_beta': 0.7,
	'mixing_ndim': 8,
	'diagonalization': 'david',
	'diago_david_ndim': 4,
	'diago_full_acc': True,
    },
    'ions': {
	'ion_dynamics':'bfgs',
    },
    'cell': {
    	'cell_dynamics':'bfgs',
    }
}

comp_params = {
	'nd':1
}

#Run Defect Tracker
from mp_api.client import MPRester as MPRester_new
with MPRester_new('API_KEY') as mp_new:
	primitive_structure = mp_new.get_structure_by_material_id(molecule_species_mpid)

"""
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
conventional_structure = SpacegroupAnalyzer(primitive_structure).get_conventional_standard_structure()
"""

tracker = Tracker.Tracker("self-defined",
                          pymatgen_unitcell = primitive_structure,
                          supercell_dim=supercell_dim,
                          molecule_species_mpid=molecule_species_mpid,
                          delete_wfc_files=True)

tracker.calculatePristineEnergy(input_data,pseudopotentials,kpts,compute_parameters=comp_params,force_recalc = force_recalc)

