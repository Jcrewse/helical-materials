from src.ase_systems import H2_chain
from gpaw import GPAW, PW, FermiDirac, hamiltonian

system = H2_chain.System(twist_angle=0, cell_size=1)

calc = GPAW(mode=PW(200),
            xc='PBE',
            kpts=(8, 1, 1),
            random=True,  # random guess (needed if many empty bands required)
            occupations=FermiDirac(0.01),
            txt='Si_gs.txt')
system.Atoms.calc = calc
system.Atoms.get_potential_energy()
ef = calc.get_fermi_level()
calc.write('Si_gs.gpw')

density = calc.density
ham = calc.hamiltonian
KE = ham.calculate_kinetic_energy(density)

# Restart from ground state and fix potential:
calc = GPAW('Si_gs.gpw').fixed_density(
    nbands=16,
    symmetry='off',
    kpts={'path': 'GX', 'npoints': 60},
    convergence={'bands': 8})

bs = calc.band_structure()
bs.plot(filename='bandstructure.png', show=True, emax=10.0)