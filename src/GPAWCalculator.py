# GPAW calculator module ######################################################
import os
import numpy as np
import matplotlib.pyplot as plt
from ase.parallel import world
from ase.units import Bohr, Ha
from ase.io import write
from gpaw import GPAW
from gpaw.unfold import Unfold, find_K_from_k, plot_spectral_function

def calc_groundstate(system, params, restart = True):
    '''
    Calculates the ground state density. 
    ------------------------------------
    Writes converged density to .gpw file.
    Calculation details written to .log file.
    '''
    
    ########## Define output filenames ########## 
    gpw_outfile = f'./{system.outname}_GS.gpw'
    log_outfile = f'./{system.outname}_GS.log'
    
    ########## Set up GPAW calculator ##########
    if os.path.isfile(gpw_outfile) and restart:
        
        if world.rank == 0:
            print(f'Ground state restart file found: {system.outname}_GS.gpw') 
        calc = GPAW(gpw_outfile)
        
        system.Atoms.calc = calc
        
    else:
        # User output of parameters
        if world.rank == 0:
            print('Calculating ground state...\n')
            print('Input parameters:')
            print(f'    calc mode: {params["mode"]}')
            print(f'    xc func:   {params["xc"]}')
            print(f'    k points:  {params["kpts"]}')
            print(f'    converge:  {params["convergence"]}')
        
        calc = GPAW(txt = log_outfile, **params)

        # Associate the GPAW calculator with the Atoms object
        system.Atoms.calc = calc

    ########### Perform some calculations ##########
    # Atomic energies
    system.e_ion_total     = system.Atoms.get_total_energy()     
    system.e_ion_potential = system.Atoms.get_potential_energy()
    system.ion_forces      = system.Atoms.get_forces()
    
    # Electronic energies
    system.e_kinetic  = calc.hamiltonian.e_kinetic
    system.e_coulomb  = calc.hamiltonian.e_coulomb
    system.e_exchange = calc.hamiltonian.e_xc
    system.e_fermi    = calc.get_fermi_level()  
    
    if not (os.path.isfile(gpw_outfile) and restart):
        calc.write(gpw_outfile, mode = 'all')
    
    ########### Output results ##########
    if world.rank == 0:
        print('\nGround state converged...\n')
        print('\n---------- Energies ----------')
        print(f'Kinetic Energy:   {Ha*system.e_kinetic:3.5f}eV')
        print(f'Potential Energy: {Ha*system.e_coulomb:3.5f}eV')
        print(f'Exchange Energy:  {Ha*system.e_exchange:3.5f}eV')
        print(f'Fermi Energy:     {system.e_fermi:3.5f}eV')
        print('\n---------- Forces ----------\n')
        print(f'Ion Forces: {system.ion_forces}')
    
    return

def calc_wavefunction(system, kpt):

    # Create folder for wavefunction files
    wf_dir = f'./{system.outname}_wfs'
    try:
        if not os.path.isdir(wf_dir): os.mkdir(wf_dir)
    except(FileExistsError):
        pass
    
    # loop over all wfs and write their cube files
    nbands = system.Atoms.calc.get_number_of_bands()
    for band in range(nbands):
        
        # Get pseudo wave functions for the band
        wf = system.Atoms.calc.get_pseudo_wave_function(band=band, kpt=kpt)
        
        # Separate real/imaginary components
        re_wf = wf.real
        im_wf = wf.imag
        
        # Calculate porbability density
        wf_2  = np.conjugate(wf)*wf
        
        # Establish file names
        re_fname  = f'./{wf_dir}/{system.outname}_{band}_RePsi.cube'
        im_fname  = f'./{wf_dir}/{system.outname}_{band}_ImPsi.cube'
        wf2_fname = f'./{wf_dir}/{system.outname}_{band}_Psi2.cube'
        
        # Write each component to files
        write(re_fname, system.Atoms,  data = re_wf)
        write(im_fname, system.Atoms,  data = im_wf)
        write(wf2_fname, system.Atoms, data = wf_2)
        
    return

def calc_bandstructure(system, npoints, unfold=False):
    '''
    Calculate the band structure of the input system
    -----------------------------
    Using fixed precalculated groundstate density read from .gpw. 
    Converged band structure written to .json file.
    Calculation details written to .log file. 
    Unfolding of supercell calculations.
    '''
    
    ########## Establish file names for calculation i/o ##########
    # gpw_infile   - GPAW restart file from density calculation
    # gpw_outfile  - GPAW restart file from bands calculation
    # png_outfile  - Bands plot output file
    # log_outfile  - Bands log output file
    # json_outfile - Bands data output file
    gpw_infile   = system.outname + '_GS.gpw'
    gpw_outfile  = system.outname + '_BS.gpw'
    png_outfile  = system.outname + '_BS.png'
    log_outfile  = system.outname + '_BS.log'
    json_outfile = system.outname + '_BS.json'

    ########## Retrieve energies for reference points #########
    e_ground = system.e_ion_potential
    e_fermi  = system.e_fermi

    ########## Retrieve unit/supercell, bandpath ##########
    cell = system.Atoms.get_cell()
    kpath = cell.bandpath(npoints = npoints, pbc = system.pbc)
    #print(len(kpath.kpts))

    ########## User output ##########
    if world.rank == 0:
        print(f'Bandpath: {kpath}')
        print([point for point in kpath.path if point != ','])

    ########### Supercell band structure unfolding ##########
    #########################################################
    if unfold:
        
        # Defines an x-axis suitable for plotting bands
        x, X, _ = kpath.get_linear_kpoint_axis()
        
        # Retrieve supercell transform matrix
        M = system.supercell_transform
        
        # Find supercell K-points from unit cell k-points
        Kpath = []
        for k in kpath.kpts:
            K = find_K_from_k(k, M)[0]
            Kpath.append(K)
            
        # Calculate band structure over supercall Brillouin zone
        bs_calc = GPAW(gpw_infile).fixed_density(
            symmetry = 'off',
            nbands = system.n_bands,
            kpts = Kpath, 
            convergence={'bands':'occupied'},
            txt = log_outfile)
        
        # Write band structure restart file
        bs_calc.write(gpw_outfile, 'all')
            
        # Initialize Unfold object
        unfolded = Unfold(name = f'{system.outname}_unfolded',
                          calc = gpw_outfile,
                          M = M,
                          spinorbit = False)
        
        special_points = [point for point in kpath.path if point != ',']
        # Calculate weights and spectral function
        # Produces two outputs:
        # 'weights_{system.outname})_unfolded.pckl'
        #       Contains \epsilon_Km and P_Km(k)
        # 'sf_{system.outname}_unfolded.pckl'
        #       Contains spectral function and coorseponding energy array
        unfolded.spectral_function(kpts = kpath.kpts, 
                                 x = x,
                                 X = X,
                                 points_name=special_points)
        
        # Plot spectral function
        plot_spectral_function(filename = f'sf_{system.outname}_unfolded',
                       eref = system.e_fermi,
                       emin = -15,
                       emax = 15)
        
    ########## Standard Band Structure Calculation ##########
    #########################################################
    else:
        # Converge the band structure non self-consistently, with a fixed density
        bs_calc = GPAW(gpw_infile).fixed_density(
            symmetry = 'off',
            nbands = system.n_bands,
            kpts = kpath, 
            convergence={'bands':'occupied'},
            txt = log_outfile)
        
        # Write band structure restart file
        bs_calc.write(gpw_outfile, 'all')

        # Create band structure object for plotting, output data
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
        # band structure axis
        bs_ax.set_title(f'System: {system.tag}, $\epsilon_0 = {e_ground:3.2f}$eV, $\epsilon_F = {e_fermi:3.2f}$eV')
        bs_ax.set_ylabel(r'$\epsilon\; [eV]$')
        bs_ax.set_ylim(system.emin, system.emax)
        
        # dos axis
        dos_ax.set_xlabel(r'$D(\epsilon) \; (total)$')
        dos_ax.set_xlim(left=0)
        dos_ax.yaxis.tick_right()
        dos_ax.yaxis.set_visible(False)
        
        plt.subplots_adjust(wspace=0.1)

        # Save figure and show
        plt.savefig(png_outfile)
        if world.rank == 0:
            plt.show()

    return

    