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
    
    def __init__(self, twist_angle = 0, cell_size = 1, layers = None, max_el = 6, sc_file = None):
        super().__init__(twist_angle, cell_size, layers)
        
        # System tag
        self.tag = 'hBN'
        
        # Lattice parameters
        self.a = self.b = 1.42
        self.c = 3.28
        
        # Supercell parameters
        self.max_el = max_el

        # Output file name
        self.outname = f'{self.tag}-Nphi-{self.N_phi}'

        # Band structure parameters
        self.k_path = 'GMKGALH'
        self.n_bands = 14
        self.emin = -7.5
        self.emax = 15
        
        self.pbc = (True, True, True)
        
        self.Atoms = self.build(self.a, self.c, self.sc_file)
        
    def build(self, a, c, sc_file):
        
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
            try:
                # Open the pickle file created from supercell.py
                pickle_name = self.sc_file
                res = pickle.load(open(pickle_name, 'rb'))
                
                # Retrieve the supercell transform matrix for supercell layer
                M = res.M()
                
                # Convert to 3x3 matrix including z-axis transform
                self.supercell_transform = np.array([[M[0][0], M[0][1], 0],
                                                     [M[1][0], M[1][1], 0], 
                                                     [0, 0, self.N_phi]])
                
                # Save the supercell as a POSCAR file
                poscar_name = f'{self.sc_outname}_SC.POSCAR'
                res.superlattice().save_POSCAR(poscar_name)
                
                # Read the POSCAR file in an Atoms object, set pbc
                hbn = io.read(poscar_name, format = 'vasp')
                hbn.set_pbc((True, True, True))
                
                return hbn
            # If the pickle file has not yet been made, warn user
            except (FileNotFoundError):
               print(f'No pickle file ({self.outname}_SC.pckl) for supercell. Run supercell.py first.')

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
            self.outname = f'{self.tag}-Cell-{cell_size}'
        else: 
            self.outname = f'{self.tag}-Nphi-{self.N_phi}'
            
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
            self.outname = f'{self.tag}-Cell-{cell_size}'
        else: 
            self.outname = f'{self.tag}-Nphi-{self.N_phi}'
            
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
            self.outname = f'{self.tag}-Cell-{cell_size}'
        else: 
            self.outname = f'{self.tag}-Nphi-{self.N_phi}'
            
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