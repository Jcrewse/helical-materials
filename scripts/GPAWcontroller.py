import os, time
from ase.units import Ha
from ase.visualize import view
from ase.parallel import world
from gpaw import PW, FermiDirac
from src.GPAWCalculator import *
import src.systems as systems
from math import pi

# Calculation parameters ######################################################
# Parameters not list are GPAW calc.default_parameters
params = {
    'mode'        : PW(500),             # Calculation mode        
    'kpts'        : (5,5,5),             # k-points sampled in periodic sys
    'random'      : True,                # Random guess of WF's in empty bands
    'xc'          : 'PBE',               # Exchange-correlation function
    'occupations' : FermiDirac(0.01),    # Occupation smearing (input # = kT)
    'convergence' : {'energy' : 0.0001}  # Convergence criteria
    }
###############################################################################

# Create System ###############################################################
system = systems.hBN(twist_angle = 2*pi/4, max_el=15)
#system.show(repeat=(1,1,1)) # 

ti = time.perf_counter()
# Ground State Calculations ###################################################
if world.rank == 0: 
    print('\n========== CALCULATING GROUND STATE ==========\n')
calc_groundstate(system, params, restart=False)
    
# Calculate Wave Functions ####################################################
if world.rank == 0: 
    print('\n========== CALCULATING WAVEFUNCTIONS ==========\n')
calc_wavefunction(system, kpt=0)

# Calculate Band Structure ####################################################
if world.rank == 0: 
    print('\n========== CALCULATING BAND STRUCTURE ==========\n')
calc_bandstructure(system, npoints=50, unfold=False)

tf = time.perf_counter()
print(f'Time: {np.round(tf-ti, 3)}s')