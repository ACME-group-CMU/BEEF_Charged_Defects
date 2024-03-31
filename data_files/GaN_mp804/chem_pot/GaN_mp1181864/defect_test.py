#Pymatgen
#from pymatgen.analysis.defects.generators import VacancyGenerator
from pymatgen.ext.matproj import MPRester

#Other
import importlib
import os
import subprocess
import shutil

#Defect Tracking
from defect_tracker import Tracker
importlib.reload(Tracker)

os.environ['ESPRESSO_PW_EXE'] = "pw.x"
os.environ['ESPRESSO_PP_EXE'] = "pp.x"
os.environ['COFFEE_DIR'] = "/trace/home/atimmins/packages/coffee/CoFFEE-main/"

# STUFF TO CHANGE
molecule_species_mpid = 'mp-1181864'
structure_path = "/trace/group/acmegroup/atimmins/testing/BEEF_Research/GaN_mp804/structure_files/mp-1181864_primitive"
kpts = [[6,6,2],[6,6,2]]
ibrav = 4

#Don't change below here!!!

force_recalc = False

# Define pseudo potentials
pseudopotentials = {'Ga': 'Ga.pbe-dn-kjpaw_psl.1.0.0.UPF',
                    'N': 'N.pbe-n-kjpaw_psl.1.0.0.UPF'}

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
        'prefix':'PbI2',
        'verbosity':'high',
        'restart_mode': 'restart',
        'tstress': True,
        'tprnfor': True,
        'max_seconds': 154800,
    },
    'system' : {
        'nat': 1,
        'ntyp': 1,
        'ecutwfc':70,
        'ecutrho':350,
        'tot_charge':0,
        'ibrav': ibrav,
        'degauss':0.01,
        'smearing':'gauss',
        'occupations':'smearing',
        'input_dft':'BEEF-vdW',
    },
    'electrons': {
        'electron_maxstep':100,
        'scf_must_converge':True,
        'conv_thr': 1e-7,
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

import json
from pymatgen.core import Lattice, Structure

with open(structure_path, 'r') as f:
    d = json.load(f)
    my_structure = Structure.from_dict(d)

tracker = Tracker.Tracker(unit_cell_type = "self-defined",
                                pymatgen_unitcell= my_structure,
                                supercell_dim = supercell_dim,
                                molecule_species_mpid=molecule_species_mpid,
                                delete_wfc_files=False)

tracker.calculatePristineEnergy(input_data,pseudopotentials,kpts,compute_parameters=comp_params,force_recalc = force_recalc)

