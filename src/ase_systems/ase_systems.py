from ase import Atoms
from ase.visualize import view
from math import pi, cos, sin

'''TO DO: INTENDED TO GENERALIZE CHAIN SYSTEMS TO A SINGLE CLASS. LOW PRIORITY.
'''

class ChainSystems():
    '''
    Class defining a number of different twisted molecular chains in ASE.
    ------------------------------------------------
    Description:
        Creates an ase.Atoms object and relevant attributes for the chain 
        of interest.
    '''

    # Constructor
    def __init__(self, twist_angle = 0, cell_size = 1):
        
        # Chain systems will all have 1D periodic boundary conditions
        self.pbc = (True, False, False)

        # System parameters
        self.twist_angle = twist_angle
        self.cell_size   = cell_size
        
        # Energies
        self.e_fermi = 0
        self.e_pot   = 0
        
        # Determine Number of layers in supercell
        if twist_angle == 0:
            N_phi = self.N_phi = cell_size
        else:
            N_phi = self.N_phi = int(2*pi/twist_angle)
        
        # Primitive cell --> Super cell transform matrix
        self.supercell_transform = [[N_phi,0,0],[0,1,0],[0,0,1]]

        self.emin = -45
        self.emax = 45
        
        
    # Hydrogen molecule chain #################################################
    def create_H2_chain(self):

        # Identifier for various usage
        self.tag = 'H2'

        # Bong lengths
        d = 0.779 # H-H bond length (determined from relaxation)
        l = 1.1   # H2-H2 bong length (value used in Hoffman chain)
            
        # Set number of bands to calculate
        # Is there a better way to determine the number of bands? 
        self.n_bands = 10*self.N_phi
            
        # Unit cell sizes
        a = self.N_phi*l
        b = 5
        c = 5

        # Create the filename string
        self.outname = f'{self.tag}_Nphi-{self.N_phi}'

        # Create a generator of angles for producing atomic positions

        # Calculate atomic positions
        positions = []
        atoms = []
        angles = [n*self.twist_angle for n in range(self.N_phi)]
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

        # Initialize Atoms object as instance variable
        self.Atoms =  Atoms(atoms,
                            positions = positions,
                            cell = [(a, 0, 0),
                                    (0, b, 0),
                                    (0, 0, c)],
                            pbc = self.pbc)
        
        return
    
    # C2H2 Chain ##############################################################
    def create_C2H2_chain(self, angle, cell_size):
        '''
        =======================================================================
        Creates an infinite chain of C2H2 molecules
        ------------------------------------------------
        Description:
            Creates and ase.Atoms object of the unit cell of an infinite chain
            of C2H2 molecules. A twist angle along the z-axis is also
            available. 

        Args:
            angle   (float): Twist angle along the x-axis
            cell_size (int): Number of unit cells to consider in calculation.
        =======================================================================
        '''
        self.tag = 'C2H2'
        # Hydrogen/Carbon position (from nist.gov data)
        h0 = 1.6644
        c0 = 0.6013
        
        # layer bong length
        l = 1.2026

        # Number of layers in the unit cell
        if angle == 0: 
            N_phi = cell_size
        else:
            N_phi = int(2*pi/angle)
        print('Number of layers in unit cell: {}'.format(N_phi))

        # Unit cell sizes
        a = N_phi*l
        b = 5
        c = 5

        # Brillouin zone path
        self.k_path = 'GX'

        # Band structure parameters
        self.n_bands = 20*N_phi
        self.emin = -25
        self.emax = 50

        self.outname = 'C2H2-pi-{}_H{:3.2f}_C{:3.2f}_l{:3.2f}_a{:3.2f}_c{:3.2f}'\
                       .format(N_phi,h0,c0,l,a,c)

        # Generator of angles for producing positions along twist axis
        angles = [n*angle for n in range(N_phi)]

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
                    pbc = (True, False, False))

    
    def show(self):
        view(self.Atoms)
        return