#Pymatgen
from pymatgen.analysis.defects.generators import VacancyGenerator
from pymatgen.ext.matproj import MPRester

#Other
import importlib
import os
import subprocess

#Defect Tracking
from qe_defect_tracker import Tracker
importlib.reload(Tracker)

os.environ['ESPRESSO_PW_EXE'] = "pw.x"
os.environ['ESPRESSO_PP_EXE'] = "pp.x"
os.environ['COFFEE_DIR'] = "/trace/home/atimmins/packages/coffee/CoFFEE_2.0/"

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
supercell_dim = [5,5,5]

# Define the input data specific to Quantum Espresso
input_data = {
    'control' : {
        'calculation': 'ensemble',
        'outdir':'./',
        'pseudo_dir': '/trace/group/acmegroup/atimmins/packages/espresso/pseudo/psl_1.0.0/pslibrary.1.0.0/pbe/PSEUDOPOTENTIALS/',
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
	'electron_maxstep':100,
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

"""
#Create a pristine supercell
tracker = Tracker.Tracker("conventional",supercell_dim=supercell_dim,molecule_species_mpid=molecule_species_mpid,delete_wfc_files=False)
tracker.calculatePristineEnergy(input_data,pseudopotentials,kpts,comp_params,force_recalc=force_recalc)
"""

from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read as ase_read

output_file = "/trace/group/acmegroup/atimmins/testing/BEEF_Research/V_C/experiments/ex30_5x5x5/supercells/PRISTINE_0/0/unitcell/espresso.pwo" 
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
		-2: 1998,
		-1: 1998,
       	    	0: 1998,
		1: 1998,
		2: 1998,
       		},
	'sigma_by_charge':{
        	},
	'epsilon':5.829713979999983, 
	'grid_size':[216,216,216],
  	},
}

tracker.calculatePristineEnergy(input_data,pseudopotentials,kpts,comp_params,correction_params=correction_params,force_recalc=force_recalc)


#Add defects
defect = {'C': -1}
relevant_species = {}
charge_list = [-2,-1,0,1,2]
tracker.addDefects("VAC", defect,charge_list)

#Calculate Defect Energies
tracker.defect_dict
tracker.getDefectEnergies()

tracker.calculateDefectEnergy(force_recalc=force_recalc,correction_params=correction_params)

"""

energy_dict = AlN_defects.getDefectEnergies()
print(energy_dict)

defect_dict = AlN_defects.getPymatgenDefectObjects()
print(defect_dict)

from pymatgen.core import Element
cp = {Element("Al"): -99.7202664,
      Element("N"): -273.5905338,
    Element("Hf"): -6088.295721,
    Element("Zr"): -1291.610034,
    Element("La"): -869.4946048,
    Element("Y"): -1050.857411}

#Plot results
#%matplotlib inline
from pymatgen.analysis.defects.thermodynamics import DefectPhaseDiagram

dentries = AlN_defects.getPymatgenDefectObjects()
dpd = DefectPhaseDiagram(dentries[0],vbm = -16.9696,band_gap = 50,filter_compatible=False)
p = dpd.plot(cp, ylim=None, ax_fontsize=1.3, fermi_level = 0)
p.show()
p.savefig("plot.png")


"""
