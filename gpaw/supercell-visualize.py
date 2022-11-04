import supercell_core as sc
from ase import io
from ase.visualize import view
from math import sqrt
import os

os.system('cls' if os.name == 'nt' else 'clear')

# hBN lattice parameters
a = 1.42
c = 6.709

# Twist angle (in degrees)
phi = 13.0

# Number of layers in structure
N_layers = 2

# User output
print('=====+++++Calculating supercell+++++=====')
print('N layers:    {:2d}'.format(N_layers))
print('Twist angle [deg]: {:1.2f}'.format(phi))
print('Twist angle [rad]: {:1.5f} '.format(phi*sc.DEGREE))

# hBN is isotructural to graphene. Create graphene like lattice.
graphene = sc.lattice()
graphene.set_vectors([3*a/2, a*sqrt(3)/2, 0], [3*a/2, -a*sqrt(3)/2, 0], [0, 0, c])
graphene.add_atom("B", (0, 0, 0)).add_atom("N", (a, 0, 0))

# Begin creating the structure
# Create a substrate (first layer)
structure = sc.heterostructure().set_substrate(graphene)

# Add layers on top of substrate
layer_angles = []
for n in range(N_layers-1):
    structure.add_layer(graphene)
    layer_angle = (n+1)*phi
    layer_angles.append([layer_angle*sc.DEGREE])

res = structure.opt(max_el = 6, thetas=layer_angles, log = True)

print('N atoms: ' + str(res.atom_count()))

res = structure.calc(M = res.M(), thetas = res.thetas())

res.superlattice().save_POSCAR('POSCAR')

atoms = io.read('POSCAR', format = 'vasp')

view(atoms, repeat = (3,3,1))

