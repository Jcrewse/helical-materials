from gpaw import GPAW, PW, FermiDirac
from gpaw.unfold import Unfold, find_K_from_k, plot_spectral_function
from src.ase_systems import H2_chain

# Initialize system
system = H2_chain.System(cell_size=2, l = 1.5)

# SCF density convergence
params = {
    'mode'        : PW(500),            # Calculation mode
    'kpts'        : (15,1,1),            # k-points sampled in periodic sys
    'random'      : True,                # Random guess of WF's in empty bands
    'xc'          : 'PBE',               # Exchange-correlation function
    'occupations' : FermiDirac(0.01),    # Occupation smearing (input # = kT)
    'convergence' : {'energy' : 0.0005}  # Convergence criteria
    }

calc = GPAW(txt='unfold_test.txt', **params)

system.Atoms.calc = calc
system.Atoms.get_potential_energy()
calc.write('unfold_test.gpw')

# Band path calculations
pc = system.Atoms.get_cell(complete = True)
bp = pc.bandpath(npoints=48, pbc = system.pbc)
x, X, _ = bp.get_linear_kpoint_axis()

M = system.supercell_transform

Kpts = []
for k in bp.kpts:
    K = find_K_from_k(k, M)[0]
    Kpts.append(K)

# Band structure calculation
calc_bands = GPAW('unfold_test.gpw').fixed_density(
    kpts=Kpts,
    symmetry='off',
    nbands=20,
    convergence={'bands': 'occupied'})

calc_bands.write('bands_unfold_test.gpw', 'all')

# Now perform the unfolding
unfold = Unfold(name='unfolded',
                calc='bands_unfold_test.gpw',
                M=M,
                spinorbit=False)

# unfold.spectral_function(kpts=bp.kpts, x=x, X=X,
#                          points_name=['G', 'X'])

# plot_spectral_function(filename = 'sf_unfolded',
#                        eref = system.e_fermi,
#                        emin = -15,
#                        emax = 15)