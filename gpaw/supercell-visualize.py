import supercell_core as sc
from ase import io
from ase.visualize import view
from math import sqrt

a = 1.42
c = 3.28
graphene = sc.lattice()
graphene.set_vectors([3*a/2, a*sqrt(3)/2, 0], [3*a/2, -a*sqrt(3)/2, 0], [0, 0, c])
graphene.add_atom("B", (0, 0, 0)).add_atom("N", (a, 0, 0))

structure = sc.heterostructure().set_substrate(graphene).add_layer(graphene)

res = structure.opt(thetas=[10*sc.DEGREE], log = True)
print(res.log)

res = structure.calc(M = res.M(), thetas = res.thetas())

res.superlattice().save_POSCAR('POSCAR')

atoms = io.read('POSCAR', format = 'vasp')

view(atoms, repeat = (3,3,1))

