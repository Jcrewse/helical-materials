import numpy as np
from gpaw import GPAW, PW
import matplotlib.pyplot as plt
from src.GPAWCalculator import calc_groundstate
from src.ase_systems import H2_chain

tolerance = 1e-3

system = H2_chain.System(twist_angle=np.pi/3)

# Convergence versus k points
pw_set = np.array([i for i in range(300, 1000, 50)])
ground_energies = []
fermi_energies = []
N_pw = []

print('\n========== Convergence Test: plane waves =========')
print('==================================================\n')
for i, n in enumerate(pw_set):
    
    params = {'mode': PW(n),
              'kpts': (50,1,1)}
    print(f'plane waves: {n}')
    
    calc_groundstate(system, params)
    
    print(f'Ground state energy: {system.e_ground:3.4f} eV')
    print(f'Fermi energy:        {system.e_fermi:3.4f} eV')
    ground_energies.append(system.e_ground)
    fermi_energies.append(system.e_fermi)
    N_pw.append(n)
    if i>0:
        delta_ground = abs(ground_energies[i] - ground_energies[i-1])
        delta_fermi  = abs(fermi_energies[i] - fermi_energies[i-1])
        if delta_ground < tolerance and delta_fermi < tolerance:
            print(f'Converged at PW# = {n}')
            break

plt.plot(N_pw, np.abs(ground_energies), label = '|Potential Energy|')
plt.plot(N_pw, np.abs(fermi_energies), label = '|Fermi energy|')
plt.xlim(N_pw[0], N_pw[-1])
plt.title(f'System: {system.tag}')
plt.xlabel('# of k points')
plt.ylabel('Energy [eV]')
plt.yscale('log')
plt.legend()
plt.savefig(f'{system.outname}_pw-convergence.png')
