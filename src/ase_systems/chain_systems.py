from math import pi, sin, cos
from ase import Atoms
from ase.visualize import view

class Chain():
    '''
    Parent class for a 1D infinite chain of molecules.
    Contains attributes present in every 1D translationally invariant system.
    '''
    
    def __init__(self, twist_angle = 0, cell_size = 1):

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
        self.pbc             = (True, False, False) # Periodic boundary conditions
        self.twist_angle     = twist_angle
        self.cell_size       = cell_size
        
        # Number of layers in the super cell
        if twist_angle == 0: 
            self.N_phi = cell_size
        else:
            self.N_phi = int(2*pi/twist_angle)
            
    def show(self):
        view(self.Atoms)
        return
        
# H2 Chain ####################################################################
class H2_chain(Chain):
    '''
    Class defining a chain of hydrogen molecules.
    Default bond lengths determined from R. Hoffman. 
    '''
    
    def __init__(self, d=0.779, l=1.5, twist_angle=0, cell_size=1):
        super().__init__(twist_angle, cell_size)
    
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
class CO_chain(Chain):
    '''
    Class defining a chain of CO molecules.
    Default bond lengths determined from relaxation calculations.
    '''
    
    def __init__(self, d=1.056, l=1.5, twist_angle=0, cell_size=1):
        super().__init__(twist_angle, cell_size)
    
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
class C2H2_chain(Chain):
    '''
    Class defining a chain of CO molecules.
    Default lengths determined from nist.gov data.
    '''
    
    def __init__(self, c0=0.779, h0=1.5, l=1.2026, twist_angle=0, cell_size=1):
        super().__init__(twist_angle, cell_size)
    
        # Energy limits for band structure plot
        self.emin    = -25
        self.emax    = 50
        self.n_bands = 15*self.N_phi
            
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
            
        self.Atoms = self.build(a, b, c, c0, h0, twist_angle)
            
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