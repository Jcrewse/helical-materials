import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from gpaw import GPAW, PW, FermiDirac
from src.ase_systems import H2_chain
from src.ase_systems import CO_chain

params = {
    'mode'        : PW(1000),            # Calculation mode
    'kpts'        : (25,1,1),            # k-points sampled in periodic sys
    'random'      : True,                # Random guess of WF's in empty bands
    'xc'          : 'PBE',               # Exchange-correlation function
    'occupations' : FermiDirac(0.01),    # Occupation smearing (input # = kT)
    'convergence' : {'energy' : 0.0001}  # Convergence criteria
    }

ls = []
es = []
for l in np.linspace(4.0,4.8,20):
    
    # Instantiate System for current bond length
    system = CO_chain.System(twist_angle = np.pi, cell_size = 1, l = l)
    #system.show()
    
    # Associate GPAW calculator with Atoms
    system.Atoms.calc = GPAW(txt=f'{system.outname}_chain-relax_{l}.log', **params)
    
    # Converge density, get potential enegy of atoms
    e_pot = system.Atoms.get_potential_energy()
    
    # Append valus to lists, output to user
    ls.append(l)
    es.append(e_pot)
    print(f'Bond length:      {l:3.4f} Ang.')
    print(f'Potential energy: {e_pot:3.4f} eV')
  
# Write results to file 
with open(f'{system.outname}_chain-relax.dat', 'w') as file:
    for l, e in zip(ls,es):
        file.write('{}  {}\n'.format(l,e))
   
# Fit data to fine grid polynomial, find minimum
fit_coef = np.polyfit(ls, es, deg=4)
poly = np.polynomial.Polynomial(fit_coef)
min_l = minimize(poly, x0 = 4)
print(min_l)
   
# Plot results
plt.xlabel('l')
plt.ylabel('E')
plt.xlim(ls[0], ls[-1]) 
plt.plot(ls, es)
plt.savefig(f'{system.outname}_chain_relax.png')
plt.show()