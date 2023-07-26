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

        def __init__(self, variable1, variable2):

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

        def build(self, build_variables):

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
When creating a 3D system (e.g. twisted hexagonal boron nitride) the system creation process becomes a bit more complicate. For a system with $N_\phi$ layers, we need to create a supercell of $N_\phi$ layers that may be used as the new unit cell of the unit cell of the twisted system. In this software, we use a previously built package known as [`supercell-core`](https://github.com/tnecio/supercell-core) to calculate the lattice vectors of the supercell. See the linked GitHub for access to the relevant paper and examples of package usage. Here I will only describe the functionality of our script `supercell.py` and how it interacts with the rest of the package.

The supercell determination process is unfortunately very memory intensive in `supercell-core` and runs into issues when we attempt to calculate the supercell at runtime on Carbon. Therefore, we have offloaded the supercell creation process from the main script and it is intended to be run locally. 