# BEEF_Charged_Defects
An investigation of the utility of Bayesian Error Estimation Functionals for Charged Point Defect Calculations

The folder structure is roughly as follows:
```
1. data_files
    1. {Material Composition}_{Material ID from The Materials Project} (e.g. 'GaN_mp804')
        1. analysis - location for chemical potential data and final processed data
        2. chem_pot - location of QE calculations for phases that are required to compute chemical potenital
        3. structure_files - location for easily referencible Pymatgen material structures
        4. {Defect Structure} (e.g. 'V_N','Ga_N','I_PB')
            1. experiments
                1. {experiment ID} (e.g. 'ex17_GaN')
                    1. defect_test.py - Experiment conditions
                    2. debug.txt - All Debug Info (In theory)
                    3. supercells
                        1. PRISTINE_0 - Pristine supercell calculations
                        2. {Defect type and index} (e.g. 'VAC_1')
                            1. {Defect Charge State}
                                1. All data files for defect charge state
2. Utilites
  1. *process_data.ipynb - Used to generate final data files for analysis
  2. Quantify_relax_within_radius.ipynb - Calculates the relaxation in defect systems
```

This repository relies upon the code in:
1. qe-defect-tracker (See https://github.com/ACME-group-CMU/qe-defect-tracker)
2. CPLAP_with_BEEF (See https://github.com/ACME-group-CMU/CPLAP_with_BEEF)