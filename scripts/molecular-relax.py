from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton
from math import sqrt, pi
from src.ase_systems import H2_chain

# Script for relaxing the atoms positions of molecules ########################

# Defining the atoms in the system
system = Atoms('H2', positions=[[0.0, 0.0, 0.0],
                               [0.0, 0.0, 1.0]])

# Associate the EMT calculator with the Atoms object
system.calc = EMT()

# QuasiNewton optimizer for the system
opt = QuasiNewton(system, trajectory='relax.emt.traj')

# Perform the optimization
# Convergence criteria: f < 0.001 eV/Angstrom
opt.run(fmax=0.001)

# Calculate the bond length
pos = system.positions
bond_vec = abs(pos[0] - pos[1])
bond_length = sqrt(bond_vec[0]**2 + bond_vec[1]**2 + bond_vec[2]**2)
print(f'Bond length: {bond_length} Angstrom')

