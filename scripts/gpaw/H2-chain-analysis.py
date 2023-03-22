from src.GPAWCalculator import *
from src.ase_systems import H2_chain
from math import pi

params = {
    'mode'        : PW(500),                 # Calculation mode
    'xc'          : 'LDA',                   # Exchange-correlation func. ID
    'h'           : 0.01,                   # Spatial grid size
    'kpts'        : (10,1,1),                 # k-points sampled in periodic sys
    'convergence' : {'energy' : 0.0005}      # Convergence criteria
    }

sys = H2_chain.System()

calc_groundstate(sys, params=params)
#calc_wavefunction(sys, params)
#calc_bandstructure(sys, npoints=100)