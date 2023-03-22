# GPAW calculator module ######################################################
import os
import numpy as np
import matplotlib.pyplot as plt
from gpaw import GPAW, PW, FermiDirac

def calc_groundstate(system, params):
    '''
    ===========================================================================
    Calculate the ground state of the system.
    -----------------------------------------
    ===========================================================================
    '''
    print('\n========== Calculating ground state ==========\n')
    
    # Define output file names from system attribute        
    gpw_outfile = f'./{system.outname}_GS.gpw'
    log_outfile = f'./{system.outname}_GS.log'
    
    # Set up GPAW calculator
    calc = GPAW(txt = log_outfile, **params)

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

def calc_wavefunction(system, params):

    # Retrieve pseudo wave functions from calculator
    # -------------------------------------------------------------------------
    system.Atoms.calc.default_parameters = params
    wave_func = system.Atoms.calc.get_pseudo_wave_function(band=1)

    # Calculate and normalize probability density
    # -------------------------------------------------------------------------
    prob = np.abs(wave_func)**2
    prob = prob/np.sum(prob)
    
    print(np.shape(prob))

    # Calculate atom positions
    x_points = np.size(wave_func[:][0][0])
    y_points = np.size(wave_func[0][:][0])
    z_points = np.size(wave_func[0][0][:])

    h = params['h']

    cell = system.Atoms.get_cell()

    x_grid = np.linspace(0, cell[0][0], num=x_points)
    y_grid = np.linspace(0, cell[1][1], num=y_points)
    z_grid = np.linspace(0, cell[2][2], num=z_points)

    positions = system.Atoms.get_positions()

    # Average over one axis
    prob_xy = np.average(prob, axis=2)
    prob_yz = np.average(prob, axis=0)
    prob_zx = np.average(prob, axis=1)
    print(np.shape(prob_xy))

    fig, (ax_xy, ax_yz, ax_zx) = plt.subplots(3,1, figsize = (5,15))
    ax_xy.set_xlabel('x')
    ax_xy.set_ylabel('y')
    ax_xy.contourf(x_grid, y_grid, prob_xy, levels = 20)
    ax_yz.set_xlabel('y')
    ax_yz.set_ylabel('z')
    ax_yz.contourf(y_grid, z_grid, prob_yz, levels = 20)
    ax_zx.set_xlabel('z')
    ax_zx.set_ylabel('x')
    ax_zx.contourf(z_grid, x_grid, prob_zx, levels = 20)
    for pos in positions:
        ax_xy.scatter(pos[0]+h, pos[1]+h)
        ax_yz.scatter(pos[1]+h, pos[2]+h)
        ax_zx.scatter(pos[2]+h, pos[0]+h)

    plt.show()

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

    print('\n========== Calculating band structure ==========\n')
    print(kpath)
    
    # Converge the band structure non self-consistently, with a fixed density
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
    dos_ax.set_xlabel(r'$D(\epsilon) \; (total)$')
    dos_ax.yaxis.tick_right()
    dos_ax.yaxis.set_visible(False)
    plt.subplots_adjust(wspace=0.1)

    # Save figure and show
    plt.savefig(png_outfile)
    plt.show()
    
    return

    