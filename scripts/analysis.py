import os
from ase.units import Ha
from ase.visualize import view
from ase.parallel import world
from gpaw import PW, FermiDirac
from src.GPAWCalculator import *
from src.ase_systems import chain_systems
from src.ase_systems import bulk_systems
from math import pi

# Calculation parameters ######################################################
# Parameters not list are GPAW calc.default_parameters
params = {
    'mode'        : PW(500),             # Calculation mode        
    'kpts'        : (5,5,5),             # k-points sampled in periodic sys
    'random'      : True,                # Random guess of WF's in empty bands
    'xc'          : 'PBE',               # Exchange-correlation function
    'occupations' : FermiDirac(0.01),    # Occupation smearing (input # = kT)
    'convergence' : {'energy' : 0.005}  # Convergence criteria
    }
###############################################################################

# Create System ###############################################################
#system = chain_systems.CO_chain(twist_angle = pi)
system = bulk_systems.hBN(twist_angle = 2*pi/2, max_el = 6)
system.show(repeat=(1,1,1))

# Ground State Calculations ###################################################
if world.rank == 0: 
    print('\n========== CALCULATING GROUND STATE ==========\n')
calc_groundstate(system, params, restart=True)
    
# Calculate Wave Functions ####################################################
if world.rank == 0: 
    print('\n========== CALCULATING WAVEFUNCTIONS ==========\n')
calc_wavefunction(system, kpt=0)

# Calculate Band Structure ####################################################
if world.rank == 0: 
    print('\n========== CALCULATING BAND STRUCTURE ==========\n')
calc_bandstructure(system, npoints=500, unfold=False)