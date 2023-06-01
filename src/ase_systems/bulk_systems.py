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
            
        # Direct control of layer number if desired
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
    
    def __init__(self, twist_angle = 0, layers = None, max_el = 6):
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
        
        self.Atoms = self.build(self.a, self.c, max_el)
        
    def build(self, a, c, max_el):
        
        print(f'========== Constructing {self.N_phi} layer supercell... ==========\n')
        
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
        opt = structure.opt(max_el = max_el, 
                            thetas = layer_angles, 
                            algorithm = 'fast',
                            log=True)
        opt.log.to_csv(f'{self.outname}_maxel-{max_el}_SC.log', index=False)
        res = structure.calc(M = opt.M(), thetas = opt.thetas())
        
        # Output to user
        print(f'Number of atoms in supercell: {res.atom_count()}')
        print(f'')
        print(f'Maximum strain: {res.max_strain():3.2f}')
        print(f'Strain tensors: ')
        for i, tensor in enumerate(res.strain_tensors()):
            print(f'\n  Layer {i}:')
            print(tensor)

        # Output superlattice as a POSCAR
        res.superlattice().save_POSCAR(f'{self.outname}_maxel-{max_el}.POSCAR')

        # Read POSCAR into Atoms object
        hbn = io.read(f'{self.outname}_maxel-{max_el}.POSCAR', format = 'vasp')

        # Set periodic boundary conditions
        hbn.set_pbc((True, True, True))

        return hbn