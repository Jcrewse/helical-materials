from ase import Atoms, io
from ase.build import bulk
from ase.visualize import view
import supercell_core as sc
import numpy as np
from math import sqrt, sin, cos

class create_system():

    def __init__(self, tag, angle=None, layers=None, stacking=None, 
                 cell_size=None):
        self.tag = tag

        if tag == 'hBN':
            self.sys = self.create_hBN(layers, stacking)
        if tag == 'helical-hBN':
            self.sys = self.create_helical_hBN(angle)
        if tag == 'Si':
            self.sys = self.create_Si()
        if tag == 'H2':
            self.sys = self.create_H2_molecule()
        if tag == 'H2-chain':
            self.sys = self.create_H2_chain(angle, cell_size)
        if tag == 'C2H2-chain':
            self.sys = self.create_C2H2_chain(angle, cell_size)

        self.e_fermi = 0
        self.e_ground = 0

    def create_hBN(self, layers=None, stacking='AA'):
        '''
        =======================================================================
        Creates the hexagonal Boron Nitride system.
        -------------------------------------------
        Description: 
            Creates an ase.Atoms object cooresponding to the h-BN system for a 
            given number of layers or bulk. Stacking arrangement may also be 
            accounted for.
            Lattice parameters and atom positons taken from: 
                Rev. Mod. Phys. 81, 1 109-162
            Graphene and h-BN are isostructural compounds.

        Args: 
            layers      (int): Number of layers, defaults to None (bulk system)
            stacking (string): Stacking order identifier 'AA' or 'AB'

        Returns: 
            ase.Atoms object
        =======================================================================           
        '''
        # Lattice parameters
        self.a = self.b = a = 1.42
        self.c = c = 3.28

        # Band structure parameters
        self.n_bands = 14
        self.emin = -7.5
        self.emax = 15
        self.k_path = 'GMKG'

        # Create output file name
        # No listing of layers if bulk calculation
        if layers == None:
            self.outname = 'hBN-bulk{}_a{}_c{}'.format(stacking, a, c)
            pbc = (True, True, True)
        else:
            self.outname = 'hBN-{}{}_a{}_c{}'.format(stacking, layers, a, c)


        # Create the positions and atoms arrays for either stacking arrangement
        positions = []
        atoms = []
        # Bulk atoms
        if layers == None:

            # Set periodic boundary conditions
            pbc = (True, True, True)

            if stacking == 'AA':
                c_layer = c
                atoms.append('BN')
                positions.append((0, 0, 0))
                positions.append((a, 0, 0))

            if stacking == 'AB':
                c_layer = 2*c
                atoms.append('BNBN')
                positions.append((0, 0, 0))
                positions.append((a, 0, 0))
                positions.append((a, 0, c/2))
                positions.append((0, 0, c/2))

        # Finite layers atoms
        for layer in range(layers):

            pbc = (True, True, False)
            
            atoms.append('B')
            atoms.append('N')
            
            if layers == 1:
                c_layer = c
            else:
                c_layer = layer*c

            if stacking == 'AA':
                positions.append((0, 0, c_layer))
                positions.append((a, 0, c_layer))

            if stacking == 'AB':
                if layer%2 == 0:
                    positions.append((0, 0, c_layer))
                    positions.append((a, 0, c_layer))
                else: 
                    positions.append((2*a, 0, c_layer))
                    positions.append((a,   0, c_layer))

        # Create the unit cell
        hbn = Atoms(atoms,
                    positions = positions,
                    cell = [(3*a/2, sqrt(3)*a/2, 0),
                            (3*a/2, -sqrt(3)*a/2, 0),
                            (0, 0, c_layer)],
                    pbc = pbc)
        
        # Center atoms and add vacuum if a layered system
        if layers is not None:
            hbn.center(vacuum=10, axis=2)

        return hbn

    def create_helical_hBN(self, angle):
        '''
        =======================================================================
        Create the hexagonal Boron Nitride system with a helical twist angle.
        ------------------------------------------
        =======================================================================           
        '''
        # Lattice constants (fixed)
        self.a = self.b = 1.42
        self.c = 6.709

        # Output file name
        self.outname = 'helical-BN_ta{:3.2f}_a{:3.2f}_b{:3.2f}'\
                       .format(angle, self.a, self.b)

        # Band structure parameters
        self.k_path = 'GMKGALH'
        self.n_bands = 14
        self.emin = -7.5
        self.emax = 15

        # Twist angle parameters
        self.twist_angle = angle
        if angle == 0:
            print('Angle must be >0 for helical system')
            return
        n_layers = int(2*np.pi/angle)
        n_layers = 2
        layer_angles = []

        a = self.a
        c = self.c

        # Create supercell
        # Define unit cell of single layer
        lattice = sc.lattice()
        lattice.set_vectors([3*a/2, a*sqrt(3)/2, 0], [3*a/2, -a*sqrt(3)/2, 0],\
                            [0, 0, c])
        lattice.add_atom("B", (0, 0, 0)).add_atom("N", (a, 0, 0))
        
        # Initialize heterostructure
        structure = sc.heterostructure().set_substrate(lattice)

        # Add layers 
        for n in range(n_layers-1):
            structure.add_layer(lattice)
            layer_angle = (n+1)*angle
            layer_angles.append([layer_angle])

        # Optimize supercell 
        res = structure.opt(max_el = 6, thetas = layer_angles, log=True)
        #res.log(self.outname + '.sc')
        res = structure.calc(M = res.M(), thetas = res.thetas())
        print('Number of atoms in supercell: {}'.format(res.atom_count()))

        # Output superlattice as a POSCAR
        res.superlattice().save_POSCAR(self.outname + '.POSCAR')

        # Read POSCAR into Atoms object
        hbn = io.read(self.outname + '_POSCAR', format = 'vasp')

        # Set periodic boundary conditions
        hbn.set_pbc((True, True, True))

        return hbn

    def create_H2_chain(self, angle, cell_size):
        '''
        =======================================================================
        Creates an infinite chain of Hydrogen molecules
        ------------------------------------------------
        Description:
            Creates and ase.Atoms object of the unit cell of an infinite chain
            of Hydrogen molecules. A twist angle along the z-axis is also
            available. 
            Lattice parameters taken from R. Hoffman

        Args:
            angle   (float): Twist angle along the x-axis
            cell_size (int): Number of unit cells to include in calculation
        =======================================================================
        '''
        # H-H bond length
        d = 0.8
        # H2-H2 bong length
        l = 1.1

        # Number of layers in the unit cell
        if angle == 0: 
            N = cell_size
        else:
            N = int((np.pi)/angle)
        print('Size of unit cell: {}'.format(N))

        # Unit cell sizes
        a = N*l
        b = 3
        c = 3

        # Brillouin zone path
        self.supercell_transform = [[N,0,0],[0,1,0],[0,0,1]]
        self.n_band_points = 200
        self.pbc = (True, False, False)

        # Band structure parameters
        self.n_bands = 2*N
        self.emin = -20
        self.emax = 45

        data_dir = './gpaw/data/'
        self.outname = data_dir + 'H2-pi-{}_d{:3.2f}_l{:3.2f}_a{:3.2f}_c{:3.2f}'\
                       .format(N,d,l,a,c)

        # Create a generator of angles for producing atomic positions
        angles = [n*angle for n in range(N)]

        positions = []
        atoms = []
        for n, phi in enumerate(angles):
            position1 = ((2*n+1)*(l/2),  (d/2)*cos(phi) + (b/2),  (d/2)*sin(phi) + (c/2))
            position2 = ((2*n+1)*(l/2), -(d/2)*cos(phi) + (b/2), -(d/2)*sin(phi) + (c/2))
            positions.append(position1)
            positions.append(position2)
            atoms.append('H')
            atoms.append('H')

        h2_twisted  = Atoms(atoms,
                            positions = positions,
                            cell = [(a, 0, 0),
                                    (0, b, 0),
                                    (0, 0, c)],
                            pbc = self.pbc)

        return h2_twisted

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
        # Hydrogen/Carbon position (from nist.gov data)
        h0 = 1.6644
        c0 = 0.6013
        # layer bong length
        l = 1.2026

        # Number of layers in the unit cell
        if angle == 0: 
            N = cell_size
        else:
            N = int((np.pi)/angle)
        print('Number of layers in unit cell: {}'.format(N))

        # Unit cell sizes
        a = N*l
        b = 5
        c = 5

        # Brillouin zone path
        self.k_path = 'GX'

        # Band structure parameters
        self.n_bands = 2*N
        self.emin = -25
        self.emax = 50

        self.outname = 'C2H2-pi-{}_H{:3.2f}_C{:3.2f}_l{:3.2f}_a{:3.2f}_c{:3.2f}'\
                       .format(N,h0,c0,l,a,c)

        # Generator of angles for producing positions along twist axis
        angles = [n*angle for n in range(N)]

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

        c2h2 = Atoms(atoms,
                    positions = positions,
                    cell = [(a, 0, 0),
                            (0, b, 0),
                            (0, 0, c)],
                    pbc = (True, False, False))

        return c2h2

    def create_H2_molecule(self):

        d = 1.8 
        l = 3
        
        num = 1

        self.outname = 'H{}_d{:3.2f}_l{:3.2f}'.format(2*num, d, l)

        H2 = Atoms('HH',
                    positions = [(0, d/4, 0),
                                 (0, 3*d/4, 0)],
                    cell = [(l, 0, 0),
                            (0, l, 0),
                            (0, 0, l)],
                    pbc = (False, False, False))

        H4 = Atoms('HHHH',
                    positions = [(0, d/4, 0),
                                 (0, 3*d/4, 0),
                                 (l, d/4, 0),
                                 (l, 3*d/4, 0)],
                    cell = [(2*l, 0, 0),
                            (0, l, 0),
                            (0, 0, l)],
                    pbc = (False, False, False))

        H2.center()

        return H2

    def create_Si(self):
        '''
        =======================================================================
        Create standard Silicon diamond structure system
        ------------------------------------------------
        Data from materialsproject.org: mp-149
        Lattice constants (Angstroms): a = b = c = 5.47
        Using ASE built-in lattice vectors and positions.
        =======================================================================
        '''
        self.a = self.b = self.c = 5.47
        self.k_path = 'GXWKGLUWL'
        self.n_bands = 14
        self.outname = 'Si_a{:3.2f}'.format(self.a)

        self.emin = -7.5
        self.emax = 15

        Si = bulk('Si', 'diamond', self.a)

        return Si

    def view(self, range):
        '''
        =======================================================================
        View the structure created using ASE built-in view function
        -----------------------------------------------------------
        =======================================================================
        '''
        
        view(self.sys, repeat = range)
        #io.write(filename=self.outname + '.xyz', images = self.sys)

        return