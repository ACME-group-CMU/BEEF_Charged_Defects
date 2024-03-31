#Pymatgen
#from pymatgen.analysis.defects.generators import VacancyGenerator
from pymatgen.ext.matproj import MPRester

from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read as ase_read

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
os.environ['COFFEE_DIR'] = "/trace/home/atimmins/packages/coffee/CoFFEE_2.0/"

#Define relevant details
molecule_species_mpid = 'mp-22893'
force_recalc = False

#K-points (can check multiple)
kpts = [[4,4,2],[1,1,1]]  #update


# Define pseudo potentials
pseudopotentials = {'Pb':'Pb.pbe-dn-kjpaw_psl.1.0.0.UPF',
		    'I':'I.pbe-dn-kjpaw_psl.1.0.0.UPF'};

#Generate the ideal dimensions of the supercell (pymatgen uses this ASE function, could incorporate more cleanly in the future
#conf = bulk(molecule_species)
#supercell_dim = find_optimal_cell_shape(conf.cell, molecule_number, 'sc')
supercell_dim = [4,4,3]

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
        'ibrav': 4,
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

"""
from mp_api.client import MPRester as MPRester_new
with MPRester_new('API_KEY') as mp_new:
	primitive_structure = mp_new.get_structure_by_material_id(molecule_species_mpid)
"""

vc_relaxed_unitcell = "/trace/group/acmegroup/atimmins/testing/BEEF_Research/PbI2_mp22893/Pb_I/ex55_443/supercells/PRISTINE_0/0/unitcell/espresso.pwo"
new_atoms_loc = ase_read(vc_relaxed_unitcell)
primitive_structure = AseAtomsAdaptor.get_structure(atoms = new_atoms_loc)

tracker = Tracker.Tracker(unit_cell_type = "self-defined",
				pymatgen_unitcell= primitive_structure,
                                supercell_dim = supercell_dim,
                                molecule_species_mpid=molecule_species_mpid,
                                delete_wfc_files=False)

correction_params = {
    'type':'FNV',
    'FNV': {
        'auto_select_bands': False,
        'electrons_at_neutral': 2368,
        'defect_level': 'shallow', #deep or shallow
        'band_by_charge' : { #specify if auto_select_bands is true,
            -1:1150,
            0:1150,
	    1:1150,
	    -2:1150,
	    2:1150,
	    3:1150,
        },
        'sigma_by_charge': {
        },
	'epsilon':17,
	'grid_size':[216,216,250],
    },
}

tracker.calculatePristineEnergy(input_data,pseudopotentials,kpts,compute_parameters=comp_params,force_recalc = force_recalc,correction_params=correction_params)

defect = {'I':-1,'Pb': 1}
relevant_species = {}
charge_list = [0,1,2,3] #add back 1,2,3
tracker.addDefects("SUB", defect,charge_list)


tracker.calculateDefectEnergy(force_recalc=force_recalc,correction_params=correction_params)