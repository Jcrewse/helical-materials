import supercell_core as sc
import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms, io
from ase.visualize import view
from math import sqrt

twist_angles = np.arange(1*sc.DEGREE, 60*sc.DEGREE, 0.1*sc.DEGREE)

# Defining unit cell of h-BN
a = 1.42
c = 6.709
graphene = sc.lattice()
graphene.set_vectors([3*a/2, a*sqrt(3)/2, 0], [3*a/2, -a*sqrt(3)/2, 0], [0, 0, c])
graphene.add_atom("B", (0, 0, 0)).add_atom("N", (a, 0, 0))


for twist_angle in twist_angles:

    # Initialize the substrate (bottom layer)
    structure = sc.heterostructure().set_substrate(graphene)
    print('Twist angle: {:3.2f}'.format(twist_angle))

    n_layers = int(2*np.pi/twist_angle)
    print('N layers: {}'.format(n_layers))

    layer_angles = []
    for n in range(n_layers):
        structure.add_layer(graphene)
        layer_angle = n*twist_angle
        layer_angles.append([layer_angle])

    print('Layer angles: {}'.format(layer_angles))

    # Optimise with theta as a free parameter
    res = structure.opt(max_el=6, thetas=layer_angles, log=True)
    res.log.to_csv('opt_{:3.2f}.log'.format(twist_angle), index=False)

# res = structure.calc(M = res.M(), thetas = res.thetas())

# fig = res.superlattice().draw()
# plt.show()

# # Save supercell to VASP POSCAR
# res.superlattice().save_POSCAR("POSCAR")

# atoms = io.read('POSCAR', format = 'vasp')

# view(atoms, repeat=(1,1,1))