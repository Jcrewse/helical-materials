import os
from ase.visualize import view
from src.GPAWCalculator import *
from src.ase_systems import H2_chain
from src.ase_systems import CO_chain
from math import pi

# Calculation parameters ######################################################
# Parameters not list are calc.default_parameters
params = {
    'mode'        : PW(1000),            # Calculation mode
    'kpts'        : (50,1,1),            # k-points sampled in periodic sys
    'random'      : True,                # Random guess of WF's in empty bands
    'xc'          : 'PBE',               # Exchange-correlation function
    'occupations' : FermiDirac(0.01),    # Occupation number smearing (input # = kT)
    'convergence' : {'energy' : 0.0001}  # Convergence criteria
    }

# Create System ###############################################################
N_phi = 1
system = H2_chain.System(twist_angle = 0, cell_size = 1)
#system.show()

# Ground State Calculations ###################################################
print('\n========== Ground State Calculation ==========\n')

# Restart previous calculation if file exists
if os.path.isfile(f'{system.outname}_GS.gpw'):
    
    print(f'Ground state restart file found: {system.outname}_GS.gpw')
    
    system.calc = GPAW(f'{system.outname}_GS.gpw')
    system.e_ground = system.calc.get_potential_energy()
    system.e_fermi  = system.calc.get_fermi_level()
    
# Converge ground state density if no file exists
else:
    
    print('Calculating ground state...')
    print('Input parameters:')
    print(f'    calc mode: {params["mode"]}')
    print(f'    xc func:   {params["xc"]}')
    print(f'    k points:  {params["kpts"]}')
    print(f'    converge:  {params["convergence"]}')
    
    
    calc_groundstate(system, params=params)
    
    print(system.e_ground)
    print(system.e_tot)
    print(system.e_fermi)
    
# Calculate Wave Functions ####################################################
# calc_wavefunction(system, params)

# Calculate Band Structure ####################################################
calc_bandstructure(system, npoints=100)