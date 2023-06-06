from math import pi, sin, cos, sqrt
from ase import io
from ase.visualize import view
from ase.parallel import world
import numpy as np
import supercell_core as sc

class Bulk():
    '''
    Parent class for bulk materials for DFT calculation
    Contains attributes present in bulk solid. 
    '''
    
    def __init__(self, twist_angle = 0, layers = None):
        
        # Elemental tag for ID
        self.tag = None

        # Ionic energies
        self.e_ion_total     = None
        self.e_ion_kinetic   = None
        self.e_ion_potential = None
        
        # Electronic energies
        self.e_total         = None
        self.e_kinetic       = None
        self.e_coulomb       = None
        self.e_exchange      = None
        self.e_fermi         = None
        
        # Temperature/forces
        self.temperature     = None
        self.ion_forces      = None
        
        # Cell geometry parameters
        self.pbc             = (True, True, True) # Periodic boundary conditions
        self.twist_angle     = twist_angle
        self.layers          = layers
        
        # Number of layers in the super cell
        if twist_angle == 0:
            self.N_phi = 1
        else:
            self.N_phi = int(2*pi/twist_angle)
            
        # Direct control of layer number if desired
        if layers != None: 
            self.N_phi = layers
            
        # Supercell lattice transform matrix
        self.supercell_transform = None
            
    def show(self, repeat=(1,1,1)):
        view(self.Atoms, repeat=repeat)
        return
    
class hBN(Bulk):
    '''
    Class defining hexagonal Boron Nitride with helical twist angle between vdW layers.
    Description: 
        Creates an ase.Atoms object cooresponding to the h-BN system for a 
        given number of layers or bulk. Stacking arrangement may also be 
        accounted for.
        Lattice parameters and atom positons taken from: 
            Rev. Mod. Phys. 81, 1 109-162
        Graphene and h-BN are isostructural compounds.
    '''
    
    def __init__(self, twist_angle = 0, layers = None, max_el = 6):
        super().__init__(twist_angle, layers)
        
        # Lattice parameters
        self.a = self.b = 1.42
        self.c = 3.28
        
        # Supercell parameters
        self.max_el = max_el
        
        # Ionic energies
        self.e_ion_total     = None
        self.e_ion_kinetic   = None
        self.e_ion_potential = None
        
        # Electronic energies
        self.e_total         = None
        self.e_kinetic       = None
        self.e_coulomb       = None
        self.e_exchange      = None
        self.e_fermi         = None
        
        # Temperature/forces
        self.temperature     = None
        self.ion_forces      = None

        # Output file name
        self.outname = f'hBN-Nphi-{self.N_phi}'
        self.sc_outname = f'{self.outname}_maxel-{self.max_el}'

        # Band structure parameters
        self.k_path = 'GMKGALH'
        self.n_bands = 14
        self.emin = -7.5
        self.emax = 15
        
        self.Atoms = self.build(self.a, self.c)
        
    def build(self, a, c):
        
        if world.rank == 0:
            print(f'========== Constructing {self.N_phi} layer supercell... ==========\n')
        
        self.tag = 'hBN'
        
        # Create supercell
        # Define unit cell of single layer
        lattice = sc.lattice()
        lattice.set_vectors([3*a/2, a*sqrt(3)/2, 0], 
                            [3*a/2, -a*sqrt(3)/2, 0],
                            [0, 0, c])
        lattice.add_atom("B", (0, 0, 0)).add_atom("N", (a, 0, 0))
        
        # Initialize heterostructure
        structure = sc.heterostructure().set_substrate(lattice)

        # Add layers
        layer_angles = []
        for n in range(self.N_phi-1):
            structure.add_layer(lattice)
            layer_angle = (n+1)*self.twist_angle
            layer_angles.append([layer_angle])
        
        # Optimize supercell 
        opt = structure.opt(max_el = self.max_el, 
                            thetas = layer_angles, 
                            algorithm = 'direct',
                            log=True)
        opt.log.to_csv(f'{self.outname}_maxel-{self.max_el}_SC.log', index=False)
        res = structure.calc(M = opt.M(), thetas = opt.thetas())
        
        # Supercell transform matrix
        M = opt.M()
        self.supercell_transform = np.array([[M[0][0], M[0][1], 0],
                                             [M[1][0], M[1][1], 0], 
                                             [0, 0, self.N_phi]])
        
        # Output to user
        if world.rank == 0:
            print(f'Number of atoms in supercell: {res.atom_count()}')
            print(f'Maximum strain: {res.max_strain():3.4f}')
            print(f'Strain tensors: ')
            for i, tensor in enumerate(res.strain_tensors()):
                print(f'\n  Layer {i}:')
                print(tensor)

        # Output superlattice as a POSCAR
        res.superlattice().save_POSCAR(f'{self.sc_outname}.POSCAR')

        # Read POSCAR into Atoms object
        hbn = io.read(f'{self.sc_outname}.POSCAR', format = 'vasp')

        # Set periodic boundary conditions
        hbn.set_pbc((True, True, True))
        
        # Output number of atoms per layer
        self.count_layer_atoms()

        return hbn
    
    def count_layer_atoms(self):
        z = []
        
        # Open POSCAR file of supercell
        with open(f'{self.sc_outname}.POSCAR') as file:
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
        if world.rank == 0:
            print(f'Total atoms: {len(z)}')
            print(f'Atoms per layer: {layer_dict}')