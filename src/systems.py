from math import pi, sin, cos, sqrt
from ase import io, Atoms
from ase.visualize import view
from ase.parallel import world
import pickle
import numpy as np
import supercell_core as sc

class HelicalSystem():
    
    def __init__(self, twist_angle = 0, cell_size = 1, layers = None):
        
        # Elemental tag
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
        self.pbc             = None
        self.twist_angle     = twist_angle
        self.cell_size       = cell_size
        self.layers          = layers
        
        # Number of layers in the super cell
        if twist_angle == 0: 
            self.N_phi = cell_size
        else:
            self.N_phi = int(2*pi/twist_angle)
            
        if layers != None:
            self.N_phi = layers
            
        # Supercell tranformation matrix
        self.supercell_transform = None
        
    # Function to show the unit cell of the system    
    def show(self, repeat = (1,1,1)):
        view(self.Atoms, repeat)
        return

    
# Hexagonal Boron Nitride #####################################################
class hBN(HelicalSystem):
    '''
    Class defining hexagonal Boron Nitride with helical twist angle between 
    vdW layers.
    Description: 
        Creates an ase.Atoms object cooresponding to the h-BN system for a 
        given number of layers or bulk. Stacking arrangement may also be 
        accounted for.
        Lattice parameters and atom positons taken from: 
            Rev. Mod. Phys. 81, 1 109-162
        Graphene and h-BN are isostructural compounds.
    '''
    
    def __init__(self, twist_angle = 0, cell_size = 1, layers = None, max_el = 6):
        super().__init__(twist_angle, cell_size, layers)
        
        # System tag
        self.tag = 'hBN'
        
        # Lattice parameters
        self.a = self.b = 1.42
        self.c = 3.28
        
        # Supercell parameters
        self.max_el = max_el

        # Output file name
        self.outname = f'hBN-Nphi-{self.N_phi}'
        self.sc_outname = f'{self.outname}_maxel-{self.max_el}'

        # Band structure parameters
        self.k_path = 'GMKGALH'
        self.n_bands = 14
        self.emin = -7.5
        self.emax = 15
        
        self.pbc = (True, True, True)
        
        self.Atoms = self.build(self.a, self.c)
        
    def build(self, a, c):
        
        if world.rank == 0:
            print(f'========== Constructing {self.N_phi} layer supercell... ==========\n')
        
        # Create supercell
        # Define unit cell of single layer
        # lattice = sc.lattice()
        # lattice.set_vectors([3*a/2, a*sqrt(3)/2, 0], 
        #                     [3*a/2, -a*sqrt(3)/2, 0],
        #                     [0, 0, c])
        # lattice.add_atom("B", (0, 0, 0)).add_atom("N", (a, 0, 0))
        
        if (self.twist_angle == 0 and self.cell_size == 1 and self.layers == None): 
            # Create the unit cell
            return Atoms(['B', 'N'],
                        positions = [(0,0,0),
                                     (a,0,0)],
                        cell = [(3*a/2, sqrt(3)*a/2, 0),
                                (3*a/2, -sqrt(3)*a/2, 0),
                                (0, 0, c)],
                        pbc = self.pbc)
            
        else:
            # try:
            res = pickle.load(open(f'{self.sc_outname}_SC.pckl', 'rb'))
            res.superlattice().save_POSCAR('{self.outname}.POSCAR')
            print(res.M())
            hbn = io.read(f'{self.sc_outname}_SC.POSCAR', format = 'vasp')
            hbn.set_pbc((True, True, True))
            return hbn
            # except (FileNotFoundError):
            #    print(f'No pickle file ({self.outname}_SC.pckl) for supercell. Run supercell-optimize.py first.')
        
    #     # Initialize heterostructure
    #     structure = sc.heterostructure().set_substrate(lattice)

    #     # Add layers
    #     layer_angles = []
    #     for n in range(self.N_phi-1):
    #         structure.add_layer(lattice)
    #         layer_angle = (n+1)*self.twist_angle
    #         layer_angles.append([layer_angle])
        
    #     # Optimize supercell 
    #     opt = structure.opt(max_el = self.max_el, 
    #                         thetas = layer_angles, 
    #                         algorithm = 'direct',
    #                         log=True)
    #     opt.log.to_csv(f'{self.outname}_maxel-{self.max_el}_SC.log')
    #     res = structure.calc(M = opt.M(), thetas = opt.thetas())
        
    #     # Supercell transform matrix
    #     M = opt.M()
    #     self.supercell_transform = np.array([[M[0][0], M[0][1], 0],
    #                                          [M[1][0], M[1][1], 0], 
    #                                          [0, 0, self.N_phi]])
        
    #     # Output superlattice as a POSCAR
    #     res.superlattice().save_POSCAR(f'{self.sc_outname}.POSCAR', 
    #                                    silent = True)

    #     # Read POSCAR into Atoms object
    #     hbn = io.read(f'{self.sc_outname}.POSCAR', format = 'vasp')
        
    #     hbn.set_pbc((True, True, True))
        
    #     # Output to user
    #     if world.rank == 0:
    #         print('\nSupercell properties: \n')
    #         print(f'    Number of atoms in supercell: {res.atom_count()}')
    #         print(f'    Maximum strain: {res.max_strain():.8f}')
    #         print(f'    Supercell vectors: ')
    #         print('         c_1 = ' + str(res.superlattice().vectors()[0]))
    #         print('         c_2 = ' + str(res.superlattice().vectors()[1]))
    #         print(f'    Atoms per layer: {self.count_layer_atoms()}')
    #         # print(f'Strain tensors: ')
    #         # for i, tensor in enumerate(res.strain_tensors()):
    #         #     print(f'\n  Layer {i}:')
    #         #     print(tensor)

    #     return hbn
    
    # def count_layer_atoms(self):
    #     z = []
        
    #     # Open POSCAR file of supercell
    #     with open(f'{self.sc_outname}.POSCAR') as file:
    #         # Skip POSCAR header
    #         for _ in range(8):
    #             next(file)
                
    #         # Read atomic positions
    #         # If z value in range, add to count
    #         for line in file:
    #             z.append(float(line.split(' ')[-1]))

    #     # Counting instances of unique z values
    #     layer_dict = {i:z.count(i) for i in z}

    #     # Output results to user
    #     return layer_dict
###############################################################################
    
# H2 Chain ####################################################################
class H2_chain(HelicalSystem):
    '''
    Class defining a chain of hydrogen molecules.
    Default bond lengths determined from R. Hoffman. 
    '''
    
    def __init__(self, d=0.779, l=1.5, twist_angle=0, cell_size=1):
        super().__init__(twist_angle, cell_size)
    
        # System tag
        self.tag = 'H2'
    
        # Energy limits for band structure plot
        self.emin    = -20
        self.emax    = 20
        self.n_bands = 25*self.N_phi
            
        # Unit cell sizes
        a = self.N_phi*l
        b = 10
        c = 10

        # Supercell tranformation matrix
        self.supercell_transform = [[self.N_phi,0,0],[0,1,0],[0,0,1]]

        # Create the filename string
        if twist_angle == 0:
            self.outname = f'H2-Cell-{cell_size}'
        else: 
            self.outname = f'H2-Nphi-{self.N_phi}'
            
        self.pbc = (True, False, False)
            
        self.Atoms = self.build(a, b, c, d, l, twist_angle)
            
        return None
    
    def build(self, a, b, c, d, l, twist_angle):

        # Create a generator of angles for producing atomic positions
        angles = [n*twist_angle for n in range(self.N_phi)]

        # Calculate atomic positions
        positions = []
        atoms = []
        for n, phi in enumerate(angles):
            position1 = ((2*n+1)*(l/2),  
                         (d/2)*cos(phi) + (b/2),  
                         (d/2)*sin(phi) + (c/2))
            position2 = ((2*n+1)*(l/2), 
                         -(d/2)*cos(phi) + (b/2), 
                         -(d/2)*sin(phi) + (c/2))
            positions.append(position1)
            positions.append(position2)
            atoms.append('H')
            atoms.append('H')

        return Atoms(atoms,
                    positions = positions,
                    cell = [(a, 0, 0),
                            (0, b, 0),
                            (0, 0, c)],
                    pbc = self.pbc)
###############################################################################
        
# CO Chain ####################################################################
class CO_chain(HelicalSystem):
    '''
    Class defining a chain of CO molecules.
    Default bond lengths determined from relaxation calculations.
    '''
    
    def __init__(self, d=1.056, l=1.5, twist_angle=0, cell_size=1):
        super().__init__(twist_angle, cell_size)
        
        # System tag
        self.tag = 'CO'
    
        # Energy limits for band structure plot
        self.emin    = -35
        self.emax    = 20
        self.n_bands = 12*self.N_phi
            
        # Unit cell sizes
        a = self.N_phi*l
        b = 10
        c = 10

        # Supercell tranformation matrix
        self.supercell_transform = [[self.N_phi,0,0],[0,1,0],[0,0,1]]

        # Create the filename string
        if twist_angle == 0:
            self.outname = f'CO-Cell-{cell_size}'
        else: 
            self.outname = f'CO-Nphi-{self.N_phi}'
            
        self.pbc = (True, False, False)
            
        self.Atoms = self.build(a, b, c, d, l, twist_angle)
            
        return None
    
    def build(self, a, b, c, d, l, twist_angle):

        # Create a generator of angles for producing atomic positions
        angles = [n*twist_angle for n in range(self.N_phi)]

        # Calculate atomic positions
        positions = []
        atoms = []
        for n, phi in enumerate(angles):
            position1 = ((2*n+1)*(l/2),  
                         (d/2)*cos(phi) + (b/2),
                         (d/2)*sin(phi) + (c/2))
            position2 = ((2*n+1)*(l/2), 
                         -(d/2)*cos(phi) + (b/2), 
                         -(d/2)*sin(phi) + (c/2))
            positions.append(position1)
            positions.append(position2)
            atoms.append('C')
            atoms.append('O')

        return Atoms(atoms,
                    positions = positions,
                    cell = [(a, 0, 0),
                            (0, b, 0),
                            (0, 0, c)],
                    pbc = self.pbc)
###############################################################################

# C2H2 Chain ##################################################################
class C2H2_chain(HelicalSystem):
    '''
    Class defining a chain of CO molecules.
    Default lengths determined from nist.gov data.
    '''
    
    def __init__(self, c0=0.779, h0=1.5, l=1.2026, twist_angle=0, cell_size=1):
        super().__init__(twist_angle, cell_size)
        
        # System tag
        self.tag = 'C2H2'
    
        # Energy limits for band structure plot
        self.emin    = -30
        self.emax    = 30
        self.n_bands = 15*self.N_phi
            
        # Unit cell sizes
        a = self.N_phi*l
        b = 10
        c = 10

        # Supercell tranformation matrix
        self.supercell_transform = [[self.N_phi,0,0],[0,1,0],[0,0,1]]
        
        self.pbc = (True, False, False)

        # Create the filename string
        if twist_angle == 0:
            self.outname = f'C2H2-Cell-{cell_size}'
        else: 
            self.outname = f'C2H2-Nphi-{self.N_phi}'
            
        self.Atoms = self.build(a, b, c, c0, h0, l, twist_angle)
            
        return None
    
    def build(self, a, b, c, c0, h0, l, twist_angle):

        # Generator of angles for producing positions along twist axis
        angles = [n*twist_angle for n in range(self.N_phi)]

        positions = []
        atoms = []
        for n, phi in enumerate(angles):
            position_c1 = ((2*n+1)*(l/2),  c0*cos(phi) + (b/2),  c0*sin(phi) + (c/2))
            position_c2 = ((2*n+1)*(l/2), -c0*cos(phi) + (b/2), -c0*sin(phi) + (c/2))
            position_h1 = ((2*n+1)*(l/2),  h0*cos(phi) + (b/2),  h0*sin(phi) + (c/2))
            position_h2 = ((2*n+1)*(l/2), -h0*cos(phi) + (b/2), -h0*sin(phi) + (c/2))
            positions.append(position_c1)
            atoms.append('C')
            positions.append(position_c2)
            atoms.append('C')
            positions.append(position_h1)
            atoms.append('H')
            positions.append(position_h2)
            atoms.append('H')

        return Atoms(atoms,
                    positions = positions,
                    cell = [(a, 0, 0),
                            (0, b, 0),
                            (0, 0, c)],
                    pbc = self.pbc)