# GPAW calculator module ######################################################
from gpaw import GPAW, PW, FermiDirac
import matplotlib.pyplot as plt

def calc_groundstate(system, params):
    '''
    ===========================================================================
    Calculate the ground state of the system.
    -----------------------------------------
    ===========================================================================
    '''
    print('========== Calculating ground state ==========')
    
    # Define output file names from system attribute        
    gpw_outfile = system.outname + '_GS.gpw'
    log_outfile = system.outname + '_GS.log'
    
    # Set up GPAW calculator
    calc = GPAW(txt = log_outfile)

    # Initialize parameters
    calc.default_parameters = params

    # Associate the GPAW calculator with the Atoms object
    system.Atoms.calc = calc

    # Perform some calculations
    E       = system.e_ground = system.Atoms.get_potential_energy()
    E_fermi = system.e_fermi  = calc.get_fermi_level()
    E_tot   = system.Atoms.get_total_energy()
    
    # Write to gpw file
    calc.write(gpw_outfile, mode = 'all')
    
    # User output
    print('Total energy:     {:3.4f} eV'.format(E_tot))
    print('Potential Energy: {:3.4f} eV'.format(E))
    print('Fermi energy:     {:3.4f} eV'.format(E_fermi))
    print(f'Output file: {log_outfile}')
    print(f'Restart file: {gpw_outfile}')
    
    return

def calc_bandstructure(system, npoints=100):
    '''
    ===========================================================================
    Calculate the band structure.
    -----------------------------
    Using fixed precalculated groundstate density read from .gpw.
    ===========================================================================
    '''
    # Establish file names for calculation i/o
    gpw_infile   = system.outname + '_GS.gpw'
    png_outfile  = system.outname + '_BS.png'
    log_outfile  = system.outname + '_BS.log'
    json_outfile = system.outname + '_BS.json'

    # Retrieve energies for reference points
    e_ground = system.e_ground
    e_fermi  = system.e_fermi

    # Retrieve cell and find suitable band path from it
    cell = system.Atoms.get_cell()
    kpath = cell.bandpath(npoints = npoints, pbc = system.pbc)

    # Converge the band structure non self-consistently, with a fixed density
    print('Calculating band structure...')
    bs_calc = GPAW(gpw_infile).fixed_density(
        symmetry = 'off',
        kpts = kpath, 
        convergence={'bands':'occupied'},
        txt = log_outfile)

    # Create band structure object for plotting
    band_struct = bs_calc.band_structure()
    band_struct.write(json_outfile)

    # Get DOS and Fermi energy
    e, dos = bs_calc.get_dos()

    # Create custom figure/axes for BS and DOS
    fig, (bs_ax, dos_ax) = plt.subplots(1,2,gridspec_kw={'width_ratios':[2,1]},
                                        figsize = (10,8), sharey = True)

    # Plot band structure onto BS axis
    band_struct.plot(ax = bs_ax)

    # Plot onto DOS axis
    dos_ax.plot(dos, e)

    # Plot configuration
    bs_ax.set_title(r'System: {}, $\epsilon_0 = {:3.2f}$eV, $\epsilon_F = {:3.2f}$eV'.format(system.tag, e_ground, e_fermi))
    bs_ax.set_ylabel(r'$\epsilon\; [eV]$')
    bs_ax.set_ylim(system.emin, system.emax)
    dos_ax.set_xlabel(r'$D(\epsilon) \; (total)$')
    #dos_ax.set_ylim(system.emin,system.emax)
    dos_ax.set_xlim(0, 1.05*max(dos))
    dos_ax.yaxis.tick_right()
    dos_ax.yaxis.set_visible(False)
    plt.subplots_adjust(wspace=0.1)

    # Save figure and show
    plt.savefig(png_outfile)
    #plt.show()
    
    return

    