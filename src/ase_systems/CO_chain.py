from ase import Atoms
from ase.visualize import view
from math import pi, cos, sin

class System():
    '''
    =======================================================================
    Class for an ASE Atoms object for a chain of carbon monoxide molecules
    ------------------------------------------------
    Description:
        Creates an ase.Atoms object of the unit cell of an infinite chain
        of carbon monoxide molecules.
        Molecular bond length calculated from relaxation.py

    Args:
        angle   (float): Twist angle along the x-axis
        cell_size (int): Number of unit cells to include in calculation
    =======================================================================
    '''
    
    def __init__(self, twist_angle, cell_size = 1):
        self.tag = 'CO'
        self.twist_angle = twist_angle
        self.cell_size = cell_size
        self.e_fermi = 0
        self.e_ground = 0
        self.pbc = (True, False, False)
        self.emin = -55
        self.emax = 45
        
        self.Atoms = self.create(twist_angle, cell_size)
    
    def create(self, angle, cell_size):

        # C-O bond length
        d = 1.056

        # CO-CO bong length
        l = 1.5

        # Number of layers in the unit cell
        if angle == 0: 
            N_phi = cell_size
        else:
            N_phi = int(2*pi/angle)
            
        self.n_bands = 8*N_phi
            
        # Unit cell sizes
        a = N_phi*l
        b = 5
        c = 5

        self.supercell_transform = [[N_phi,0,0],[0,1,0],[0,0,1]]

        # Create the filename string
        self.outname = f'CO-Nphi-{N_phi}'

        # Create a generator of angles for producing atomic positions
        angles = [n*angle for n in range(N_phi)]

        # Calculate atomic positions
        positions = []
        atoms = []
        for n, phi in enumerate(angles):
            position1 = ((2*n+1)*(l/2),  (d/2)*cos(phi) + (b/2),  (d/2)*sin(phi) + (c/2))
            position2 = ((2*n+1)*(l/2), -(d/2)*cos(phi) + (b/2), -(d/2)*sin(phi) + (c/2))
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
    
    def show(self):
        view(self.Atoms)
        return