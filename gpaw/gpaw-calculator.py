'''
###############################################################################
GPAW calculator
###############################################################################
'''
import systems
import matplotlib.pyplot as plt
import os
import numpy as np
from gpaw import GPAW, PW, FermiDirac

def main():

    # Default parameters for the ground state calculation
    params = {
    'xc'          : 'PBE',                   # Exchange-correlation func. ID
    'kpoints'     : (8,1,1),                 # k-points sampled in periodic sys
    'random'      : True,                    # Random guess (of what?)
    'occupations' : FermiDirac(0.01),        # Occupation number smearing
    'convergence' : {'energy' : 0.0005}      # Convergence criteria
    }

    # Create the system of interest
    system = systems.create_system('H2-chain', angle=0, cell_size = 2)

    # User output of system for inspection
    system.view(range = (1,1,1))

    # Calculate the electronic ground state of the system
    if not os.path.isfile(system.outname + '_GS.gpw'):
        calc_groundstate(system, params = params)
    else:
        print('Restarting from: ' + system.outname)
        system.sys.calc = GPAW(system.outname + '_GS.gpw')
        system.e_ground = system.sys.calc.get_potential_energy()
        system.e_fermi  = system.sys.calc.get_fermi_level()


    # Calculate the groundstate wavefunction
    #calc_wavefunction(system, params)

    # Calculate the band structure of the system, output
    calc_bandstructure(system)

    return

def calc_groundstate(system, params):
    '''
    ===========================================================================
    Calculate the ground state of the system.
    -----------------------------------------
    ===========================================================================
    '''       
    gpw_outfile = system.outname + '_GS.gpw'
    log_outfile = system.outname + '_GS.log'

    print('Calculating ground state...')
    
    # Set up GPAW calculator
    calc = GPAW(txt = log_outfile)
    #calc = GPAW(mode = 'lcao', basis = 'dzp', txt = log_outfile)
    #calc = GPAW(mode = PW(1000), txt = log_outfile)

    calc.default_parameters = params

    # Associate the GPAW calculator with the hBN system
    system.sys.calc = calc

    # Perform some calculations
    E       = system.e_ground = system.sys.get_potential_energy()
    E_fermi = system.e_fermi  = calc.get_fermi_level()

    print('Ground state energy: {:3.4f} eV'.format(E))
    print('Fermi energy:        {:3.4f} eV'.format(E_fermi))

    # Write to gpw file
    calc.write(gpw_outfile, mode = 'all')

    return

def calc_wavefunction(system, params):

    # Retrieve pseudo wave functions from calculator
    # -------------------------------------------------------------------------
    wave_func = system.sys.calc.get_pseudo_wave_function(band=0)

    # Calculate and normalize probability density
    # -------------------------------------------------------------------------
    prob = np.abs(wave_func)**2
    prob = prob/np.sum(prob)

    # Calculate atom positions
    x_points = np.size(wave_func[:][0][0])
    y_points = np.size(wave_func[0][:][0])
    z_points = np.size(wave_func[0][0][:])

    h = params['h']

    cell = system.sys.get_cell()

    x_grid = np.linspace(0, cell[0][0], num=x_points)
    y_grid = np.linspace(0, cell[1][1], num=y_points)
    z_grid = np.linspace(0, cell[2][2], num=z_points)

    positions = system.sys.get_positions()

    z_mid = int(z_points/2)
    y_mid = int(y_points/2)
    x_mid = int(x_points/2)

    fig, (ax_xy, ax_yz, ax_zx) = plt.subplots(3,1, figsize = (8,10))
    ax_xy.set_xlabel('x')
    ax_xy.set_ylabel('y')
    ax_xy.contourf(x_grid, y_grid, prob[:][:][z_mid], levels = 20)
    ax_yz.set_xlabel('y')
    ax_yz.set_ylabel('z')
    ax_yz.contourf(y_grid, z_grid, prob[x_mid][:][:], levels = 20)
    ax_zx.set_xlabel('z')
    ax_zx.set_ylabel('x')
    ax_zx.contourf(z_grid, x_grid, prob[:][y_mid][:], levels = 20)
    for pos in positions:
        ax_xy.scatter(pos[0]+h, pos[1]+h)
        ax_yz.scatter(pos[1]+h, pos[2]+h)
        ax_zx.scatter(pos[2]+h, pos[0]+h)

    plt.show()

    return

def calc_bandstructure(system):
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
    cell = system.sys.get_cell()
    kpath = cell.bandpath(npoints = system.n_band_points, pbc = system.pbc)

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
    e, dos = bs_calc.get_dos(spin = 0, npts = 501, width = 0.1)

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
    plt.show()


    return
        
if __name__ == '__main__':
    main()