#Pymatgen
from pymatgen.analysis.defects.generators import VacancyGenerator
from pymatgen.ext.matproj import MPRester

#Other
import importlib
import os
import subprocess

#Defect Tracking
from defect_tracker import Tracker
importlib.reload(Tracker)


os.environ['ESPRESSO_PW_EXE'] = "/Users/andrewtimmins/Documents/CMU/Research/QE_ASE/dsProject/qe/qe-7.1/PW/src/pw.x"
os.environ['ESPRESSO_PP_EXE'] = "/Users/andrewtimmins/Documents/CMU/Research/QE_ASE/dsProject/qe/qe-7.1/PP/src/pp.x"
os.environ['COFFEE_DIR'] = "/Users/andrewtimmins/Documents/CMU/Research/QE_ASE/dsProject/coffee/CoFFEE_2.0"

#Variables to change


#Define relevant details
force_recalc = False
molecule_species_mpid = 'mp-66'


#K-points (can check multiple)
kpts = [[4,4,4],[1,1,1]] #update


# Define pseudo potentials
pseudopotentials = {'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF'}

#Generate the ideal dimensions of the supercell (pymatgen uses this ASE function, could incorporate more cleanly in the future
#conf = bulk(molecule_species)
#supercell_dim = find_optimal_cell_shape(conf.cell, molecule_number, 'sc')
supercell_dim = [4,4,4]

# Define the input data specific to Quantum Espresso
input_data = {
    'control' : {
        'calculation': 'ensemble',
        'outdir':'./',
        'pseudo_dir': '//group///packages/espresso/pseudo/psl_1.0.0/pslibrary.1.0.0/pbe/PSEUDOPOTENTIALS/',
        'prefix':'diamond',
        'verbosity':'high',
        'restart_mode': 'restart',
	'tstress': True,
	'tprnfor': True,
    },
    'system' : {
        'nat': 1,
        'ntyp': 1,
        'ecutwfc':70,
	'ecutrho':350,
        'tot_charge':0,
        'ibrav': 1,
	'degauss':0.01,
	'smearing':'gauss',
	'occupations':'smearing',
	'input_dft': 'BEEF-vdW',
    },
    'electrons': {
	'electron_maxstep':80,
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

comp_params = {'nd':1}

#Create a pristine supercell
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read as ase_read

output_file = "/Users/andrewtimmins/Documents/CMU/Research/QE_ASE/dsProject/Testing/BEEF_Research/V_C/experiments/ex31_4x4x4/supercells/PRISTINE_0/0/unitcell/espresso.pwo"
new_atoms_loc = ase_read(output_file) 
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
       	'electrons_at_neutral': 1022,
       	'defect_level': 'shallow', #deep or shallow
       	'band_by_charge' : {
		-2: 1022, #check the pwo's directly to get 1022
		-1: 1022,
       	    	0: 1022,
		1: 1022,
		2: 1022,
       		},
	'sigma_by_charge':{
        	},
	'epsilon':5.829713979999983,
	'grid_size':[180,180,180],
   	},
}
tracker.calculatePristineEnergy(input_data,pseudopotentials,kpts,comp_params,correction_params=correction_params,force_recalc=force_recalc)


#Add defects
defect = {'C': -1}
relevant_species = {}
charge_list = [0,-2,2,-1,1]
tracker.addDefects("VAC", defect,charge_list)

#Calculate Defect Energies
tracker.defect_dict
tracker.getDefectEnergies()

tracker.calculateDefectEnergy(force_recalc=force_recalc,correction_params=correction_params)


#%%
