# Helical Materials
Investigating nanoscale materials with helical symmetry via density functional theory simulation. Scripts to construct and calculate band structure properties of helical materials via the [GPAW](https://wiki.fysik.dtu.dk/gpaw/) density functional theory package. 

Note: This package is designed for use on the [Argonne Carbon Cluster](https://wiki.anl.gov/cnm/HPC/Carbon_Cluster_-_Overview).

## Installation
1. Close repository: `git clone https://github.com/jcrewse/helical-materials`
2. Install dependencies:  
GPAW `conda install -c conda-forge gpaw` (the pre-installed GPAW on Carbon is not compatible)  
supercell-core `pip3 install supercell-core --user`
3. Install package (from): `pip3 install -e .`

## Usage
The main script is `controller.py` from which you may select the system you would like to simulate, set the relevant parameters, and choose the calculations you would like to perform. The list of parameters in `controller.py` is not exhaustive. For a full description see: [GPAW basics](https://wiki.fysik.dtu.dk/gpaw/documentation/basic.html)

## Internals
The main portion of the calculations are contained in `GPAWcalculator.py`. This file contains functions that when called from `controller.py`, perform the relevant steps in GPAW to carry out the calculations. 

## Systems
GPAW performs calculations on ASE.Atoms objects. We define a system in the `systems.py` file. Each system is a superclass of an Atoms object with additional attributes that characterize the helical symmetries (e.g. twist angle, supercell transforms matrix, etc.) as well as relevant energies, forces, and system identifies tags.

### Example: Constructing a 1D system
Here is a template for creating a new system object in the `systems.py` module

For a 1D chain system:
```python
class YourSystem(HelicalSystem):

        def __init__(self, init_vars):

            # Label for output files, plots
            self.tag = 'YourTag'

            # Lattice constants
            self.a = #lattice_constant1
            self.b = #lattice_constant2
            self.c = #lattice_constant3

            # Bandstructure parameters
            self.
            self.emin = # Min for band structure plots
            self.emax = # Max

            # Periodic boundary conditions for 1D 
            self.pbc = (True, False, False)

            # Create the ASE.Atoms object
            self.Atoms = self.build(build_variables)

        def build(self, build_vars):

            # Create a generator of angles for producing atomic positions
            angles = [n*twist_angle for n in range(self.N_phi)]

            # Calculate atomic positions
            positions = []
            atoms = []
            for n, phi in enumerate(angles):
                position = # Positions of atoms
                positions.append(position)
                atoms.append('H') # ASE identifies atoms with their elemental symbol as a string

            return Atoms(atoms,
                positions = positions,
                cell = [(a, 0, 0),
                        (0, b, 0),
                        (0, 0, c)],
                pbc = self.pbc)
```

### Bulk systems: System creation with `supercell.py`
When creating a 3D system (e.g. twisted hexagonal boron nitride) the system creation process becomes a bit more complicate. For a system with $N_\phi$ layers, we need to create a supercell of $N_\phi$ layers that may be used as the new unit cell of the unit cell of the twisted system. In this software, we use a previously built package known as `supercell-core` to calculate the lattice vectors of the supercell. See the [GitHub](https://github.com/tnecio/supercell-core) for access to the relevant paper and examples of package usage. Here I will only describe the functionality of our script `supercell.py` and how it interacts with the rest of the package.

The supercell determination process is unfortunately very memory intensive in `supercell-core` and runs into issues when we attempt to calculate the supercell in a parallel calculation on Carbon. Therefore, we have offloaded the supercell creation process from the main script, and it is intended to be run locally prior to GPAW calculations. 

`supercell.py` determines the supercell and outputs the results in the form of three files:  

- POSCAR file: Contains atom types, positions and supercell vector information. For more information on the POSCAR formatting see the [VASP wiki](https://www.vasp.at/wiki/index.php/POSCAR).
- Pickle (.pckl) file: Serialized `heterostructure` object. Contains all information on the results of the supercell calculation. 
- Log file: Human readable output of supercell results. 

Once you have run `supercell.py` for the system of interest use the filename of the output files as an input argument to the bulk system class 

```python
system = systems.YourSystem('YourSystemFile')
```

### Example: twisted hexagonal boron nitride

```python
class hBN(HelicalSystem):
    '''
    Class defining hexagonal Boron Nitride with helical twist angle between layers.
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
        self.outname = f'{self.tag}-Nphi-{self.N_phi}'
        self.sc_outname = f'{self.outname}_maxel-{self.max_el}'

        # Band structure parameters
        self.k_path = 'GMKGALH' # Band paths as a string of relevant high-symmetry points
        self.n_bands = 14       # Number of bands to calculate
        self.emin = -7.5        # Min energy for band struct plots
        self.emax = 15          # Max
        
        # Periodic boundary conditions 
        self.pbc = (True, True, True)
        
        # Create the ASE.Atoms object
        self.Atoms = self.build(self.a, self.c)
        
    def build(self, a, c):
        
        if (self.twist_angle == 0 and self.cell_size == 1 and self.layers == None): 
            # Create the unit cell of untwisted system
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
                pickle_name = f'{self.sc_outname}_SC.pckl'
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
```

