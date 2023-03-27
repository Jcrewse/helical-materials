import numpy as np
from gpaw import GPAW, PW
import matplotlib.pyplot as plt
from src.GPAWCalculator import calc_groundstate
from src.ase_systems import H2_chain
from src.ase_systems import CO_chain

tolerance = 5e-3

system = CO_chain.System(twist_angle=2*np.pi/2)

# Convergence versus k points
Nk_set = np.array([10*(i+1) for i in range(50)])
ground_energies = []
fermi_energies = []
N_kpts = []

print('\n========== Convergence Test: k-points =========')
print('===============================================\n')
for i, n in enumerate(Nk_set):
    
    params = {'mode' : 'lcao',
              'gpts' : (n,n,n)}
    print(f'g-points: {n}')
    
    calc_groundstate(system, params)
    
    print(f'Ground state energy: {system.e_ground:3.4f} eV')
    print(f'Fermi energy:        {system.e_fermi:3.4f} eV')
    ground_energies.append(system.e_ground)
    fermi_energies.append(system.e_fermi)
    N_kpts.append(n)
    if i > 0:
        delta_ground = abs(ground_energies[i] - ground_energies[i-1])
        delta_fermi  = abs(fermi_energies[i] - fermi_energies[i-1])
        if delta_ground < tolerance and delta_fermi < tolerance:
            print(f'Converged at k-points = {n}')
            break

plt.plot(N_kpts, np.abs(ground_energies), label = 'Potential Energy')
plt.plot(N_kpts, np.abs(fermi_energies), label = 'Fermi energy')
plt.xlim(N_kpts[0], N_kpts[-1])
plt.title(f'System: {system.tag}')
plt.xlabel('# of k points')
plt.ylabel('Energy [eV]')
plt.yscale('log')
plt.legend()
plt.savefig(f'{system.outname}_k-convergence.png')
