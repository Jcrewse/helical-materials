import os
from ase.units import Ha
from gpaw import PW, FermiDirac
from ase.visualize import view
from src.GPAWCalculator import *
from src.ase_systems import H2_chain, CO_chain, C2H2_chain
from math import pi

# Calculation parameters ######################################################
# Parameters not list are GPAW calc.default_parameters
params = {
    'mode'        : PW(300),             # Calculation mode        
    'kpts'        : (30,1,1),            # k-points sampled in periodic sys
    'random'      : True,                # Random guess of WF's in empty bands
    'xc'          : 'PBE',               # Exchange-correlation function
    'occupations' : FermiDirac(0.01),    # Occupation smearing (input # = kT)
    'convergence' : {'energy' : 0.0005}  # Convergence criteria
    }
###############################################################################

# Create System ###############################################################
system = H2_chain.System(twist_angle = pi)
#system.show()

# Ground State Calculations ###################################################
print('\n========== CALCULATING GROUND STATE ==========\n')
calc_groundstate(system, params, restart=False)
    
# Calculate Wave Functions ####################################################
print('\n========== CALCULATING WAVEFUNCTIONS ==========\n')
calc_wavefunction(system, kpt=0)

# Calculate Band Structure ####################################################
print('\n========== CALCULATING BAND STRUCTURE ==========\n')
calc_bandstructure(system, npoints=50, unfold=True)