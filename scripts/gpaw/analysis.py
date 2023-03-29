import os
from ase.units import Ha
from gpaw import PW, FermiDirac
from ase.visualize import view
from src.GPAWCalculator import *
from src.ase_systems import H2_chain
from src.ase_systems import CO_chain
from src.ase_systems import C2H2_chain
from math import pi

restart = False

# Calculation parameters ######################################################
# Parameters not list are GPAW calc.default_parameters
params = {
    'mode'        : PW(1000),            # Calculation mode
    'kpts'        : (30,1,1),            # k-points sampled in periodic sys
    'random'      : True,                # Random guess of WF's in empty bands
    'xc'          : 'PBE',               # Exchange-correlation function
    'occupations' : FermiDirac(0.01),    # Occupation smearing (input # = kT)
    'convergence' : {'energy' : 0.0001}  # Convergence criteria
    }
###############################################################################

# Create System ###############################################################
N_phi = 1
system = H2_chain.System(twist_angle = 0, cell_size = 1)
#system.show()

# Ground State Calculations ###################################################
print('\n========== Ground State Calculation ==========\n')

# Restart previous calculation if file exists
if os.path.isfile(f'{system.outname}_GS.gpw') and restart:
    
    print(f'Ground state restart file found: {system.outname}_GS.gpw')
    
    system.calc = GPAW(f'{system.outname}_GS.gpw')
    system.e_pot = system.calc.get_potential_energy()
    system.e_fermi  = system.calc.get_fermi_level()
    
# Converge ground state density if no file exists
else:
    
    # User output of parameters
    print('Calculating ground state...\n')
    print('Input parameters:')
    print(f'    calc mode: {params["mode"]}')
    print(f'    xc func:   {params["xc"]}')
    print(f'    k points:  {params["kpts"]}')
    print(f'    converge:  {params["convergence"]}')
    
    # Call density convergence function
    calc_groundstate(system, params=params)
    
    # Output resulting energies
    print('\nGround state converged...\n')
    print(f'Kinetic Energy:   {Ha*system.e_kinetic:3.5f}eV')
    print(f'Potential Energy: {Ha*system.e_coulomb:3.5f}eV')
    print(f'Exchange Energy:  {Ha*system.e_exchange:3.5f}eV')
    print(f'Fermi Energy:     {system.e_fermi:3.5f}eV')
    print(f'Ion Forces: {system.ion_forces}')
    
# Calculate Wave Functions ####################################################
# calc_wavefunction(system, params)

# Calculate Band Structure ####################################################
#calc_bandstructure(system, npoints=500)