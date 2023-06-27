import supercell_core as sc
import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms, io
from ase.visualize import view
from math import sqrt
import pickle

def count_layer_atoms(poscar_file):
    z = []
    
    # Open POSCAR file of supercell
    with open(poscar_file) as file:
        # Skip POSCAR header
        for _ in range(8):
            next(file)
            
        # Read atomic positions
        # If z value in range, add to count
        for line in file:
            z.append(float(line.split(' ')[-1]))

    # Counting instances of unique z values
    layer_dict = {i:z.count(i) for i in z}

    # Output results to user
    return layer_dict
 
# User inputs
tag = 'hBN'
N_phi  = 3
max_el = 25

# Defining some lattice constants
a = 1.42
c = 3.28

# Filenames
file_name = f'{tag}-Nphi-{N_phi}_maxel-{max_el}_SC'
poscar_file = f'{file_name}.POSCAR'
log_file = f'{file_name}.log'
pickle_file = f'{file_name}.pckl'

# Defining the lattice unit cell 
# THIS NEEDS TO BE GENERALIZED TO THINGS NO ISOSTRUCTURAL TO GRAPHENE
lattice = sc.lattice()
lattice.set_vectors([3*a/2, a*sqrt(3)/2, 0], [3*a/2, -a*sqrt(3)/2, 0], [0, 0, c])
lattice.add_atom("B", (0, 0, 0)).add_atom("N", (a, 0, 0))

# Initialize the substrate (bottom layer)
structure = sc.heterostructure().set_substrate(lattice)

# Build heterostructure up layer by layer
layer_angles = []
twist_angle  = 2*np.pi/N_phi
for n in range(N_phi):
    structure.add_layer(lattice)
    layer_angle = n*twist_angle
    layer_angles.append([layer_angle])

# Optimize lattice vectors 
res = structure.opt(max_el=max_el, thetas=layer_angles, log=True)

# Save supercell to VASP POSCAR
res.superlattice().save_POSCAR(poscar_file)

# Output to user
print('\nSupercell properties: \n')
print(f'    Number of atoms in supercell: {res.atom_count()}')
print(f'    Maximum strain: {res.max_strain():.8f}')
print(f'    Supercell vectors: ')
print('         c_1 = ' + str(res.superlattice().vectors()[0]))
print('         c_2 = ' + str(res.superlattice().vectors()[1]))
print(f'    Atoms per layer: {count_layer_atoms(poscar_file)}')
print(f'    Transform matrix: M = {res.M()[0,:]}\n')
print(f'                          {res.M()[1,:]}\n')
print(f'Strain tensors: ')

# Write to log file
logfile = open(log_file, 'w')
logfile.write(f'{tag}_N_phi_{N_phi}\n')
logfile.write('Supercell properties: \n')
logfile.write(f'    Number of atoms: {res.atom_count()}\n')
logfile.write(f'    Maximum strain:  {res.max_strain():.8f}\n')
logfile.write(f'    Supercell vectors: \n')
logfile.write('         c_1 = ' + str(res.superlattice().vectors()[0]) + '\n')
logfile.write('         c_2 = ' + str(res.superlattice().vectors()[1]) + '\n')
logfile.write(f'    Atoms per layer: {count_layer_atoms(poscar_file)}\n')
logfile.write(f'    Transform matrix: M = {res.M()[0,:]}\n')
logfile.write(f'                          {res.M()[1,:]}\n')
logfile.write(f'Strain tensors: \n')

for i, tensor in enumerate(res.strain_tensors()):
    print(f'\nLayer {i}:')
    print(tensor)
    logfile.write(f'\nLayer {i}:\n')
    logfile.write(str(tensor))

logfile.close()

pickle.dump(res, open(pickle_file, 'wb'))

# View the resulting supercell
#atoms = io.read(poscar_file, format = 'vasp')
#view(atoms, repeat=(1,1,1))