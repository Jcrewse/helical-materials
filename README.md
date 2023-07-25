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
