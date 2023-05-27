from math import pi, sin, cos, sqrt
from ase import io
from ase.visualize import view
import supercell_core as sc

class Bulk():
    '''
    Parent class for bulk materials for DFT calculation
    Contains attributes present in bulk solid. 
    '''
    
    def __init__(self, twist_angle = 0, layers = None):

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
            
        if layers != None: 
            self.N_phi = layers
            
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
    
    def __init__(self, twist_angle = 0, layers = None):
        super().__init__(twist_angle, layers)
        
        # Lattice parameters
        self.a = self.b = 1.42
        self.c = 3.28

        # Output file name
        self.outname = f'hBN-Nphi-{self.N_phi}'

        # Band structure parameters
        self.k_path = 'GMKGALH'
        self.n_bands = 14
        self.emin = -7.5
        self.emax = 15
        
        self.Atoms = self.build(self.a, self.c)
        
    def build(self, a, c):
        
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
        res = structure.opt(max_el = 5, 
                            thetas = layer_angles, 
                            algorithm = 'fast',
                            log=True)
        res = structure.calc(M = res.M(), thetas = res.thetas())
        print('Number of atoms in supercell: {}'.format(res.atom_count()))

        # Output superlattice as a POSCAR
        res.superlattice().save_POSCAR(self.outname + '.POSCAR')

        # Read POSCAR into Atoms object
        hbn = io.read(self.outname + '.POSCAR', format = 'vasp')

        # Set periodic boundary conditions
        hbn.set_pbc((True, True, True))

        return hbn