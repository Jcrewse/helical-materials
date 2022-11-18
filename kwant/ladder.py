import kwant
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import scipy.sparse.linalg as sla
from math import pi, cos


hop_dict = {'t':0, 'v':1, 'w':2}

def main():

    os.system('clear')

    #analyze_bandstructure_varyt()
    #analyze_bandstructure_varyv()
    #analyze_bandstructure_varyw()
    #analyze_bandstructure_varyphi()

    # Create a finite system to visualize the structure
    system = make_system(a = 1, b = 1, phi = 0, 
                         u = 5, hops = (5,1,1), 
                         width = 2, cell_scale = 2).finalized()

    plot_band_state(system)
    #plot_spec_prob(system, n_eigs=1, plot_eigs=[0])

    return

def make_system(a, b, phi, u, hops, width, length=None, cell_scale = 1):

    t, v, w = hops

    lat = kwant.lattice.general([(a,0), (0,b)], norbs=1)

    sys = kwant.Builder()

    if length is not None:
        # Create a finite sizes system to show structure
        sys[lat.shape((lambda pos: (pos[1] >= 0 and pos[1] < b*width) and (pos[0] >= 0 and pos[0] < cell_scale*a*length)), (0,0))] = u
    else:
        # Create an infinite system for band structure calculations
        sys = kwant.Builder(kwant.TranslationalSymmetry((cell_scale*a,0)))
        sys[lat.shape((lambda pos: pos[1] >= 0 and pos[1] < b*width), (0,0))] = u

    # Intra-layer hoppings
    sys[kwant.builder.HoppingKind((0,b), lat, lat)]  = -t

    # Inter-layer hoppings
    sys[kwant.builder.HoppingKind((a,0), lat, lat)]  = -v*((0.5*a**2)*(1-cos(phi)) + b**2)**(-1)

    # Cross hoppings
    sys[kwant.builder.HoppingKind((a,b), lat, lat)]  = -w*((0.5*a**2)*(1+cos(phi)) + b**2)**(-1)
    sys[kwant.builder.HoppingKind((-a,b), lat, lat)] = -w*((0.5*a**2)*(1+cos(phi)) + b**2)**(-1)

    return sys

def plot_structure(sys):

    fig, ax = plt.subplots()

    kwant.plot(sys, ax=ax)

    ax.axis('scaled')
    ax.set_xlabel(r'$x\:[a_{}]$'.format(0))
    ax.set_ylabel(r'$y\:[a_{}]$'.format(0))
    ax.set_title('System structure')

    plt.show()

    return

def plot_bandstructure(systems, phi, hoppings, title, vary_var):

    colors = iter(cm.rainbow(np.linspace(0,1,len(systems))))

    fig, ax = plt.subplots(figsize=(6,8))
    for i, sys in enumerate(systems):
        bands = kwant.physics.Bands(sys)
        momenta = np.linspace(0, pi, 1000)
        energies = [bands(k)*CONVERT_EV for k in momenta]
        if vary_var == 'phi': label = r'$\phi = $' + str(round(phi[i],2))
        else: label = str(vary_var) + ' = {}eV'.format(round(hoppings[i][hop_dict[vary_var]], 3))
        hop_color = next(colors)
        ax.plot(momenta, energies, color=hop_color, label=label)

    ax.set_title(title)
    ax.set_xlabel(r'$k \quad [a^{-1}]$')
    ax.set_ylabel(r'$\epsilon \quad [eV]$')
    ax.set_xlim(0,pi)

    # Gets rid of duplicate labels in the legend (shamelessly copypasta'd from stackexchange)
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique), bbox_to_anchor = (1.0,1.0))

    plt.tight_layout()
    plt.show()

    return

def analyze_bandstructure_varyt():

    a      = 1
    b      = 1
    width  = 2
    phi    = 0.0
    u      = 1.0
    t_hops = [2.0, 1.0, 0.5, 0.25]
    v_hops = [1.0]
    w_hops = [0.0]
    
    hoppings = list((t,v,w) for t in t_hops for v in v_hops for w in w_hops)

    systems = []
    for hopping in hoppings:
        systems.append(make_system(a, b, phi, u, hopping, width).finalized())

    vhop_label = round(hoppings[0][1], 3)
    whop_label = round(hoppings[0][2], 3)
    title = 'Simple ladder band structure - v = {}eV, w = {}eV'.format(vhop_label, whop_label)
    plot_bandstructure(systems, phi, hoppings, title, 't')

    return

def analyze_bandstructure_varyv():

    a      = 1
    b      = 1
    width  = 2
    phi    = 0.0
    u      = 1.0
    t_hops = [1.0]
    v_hops = [2.0, 1.0, 0.5, 0.25]
    w_hops = [0.0]
    
    hoppings = list((t,v,w) for t in t_hops for v in v_hops for w in w_hops)

    systems = []
    for hopping in hoppings:
        systems.append(make_system(a, b, phi, u, hopping, width).finalized())

    thop_label = round(hoppings[0][0], 3)
    whop_label = round(hoppings[0][2], 3)
    title = 'Simple ladder band structure - t = {}eV, w = {}eV'.format(thop_label, whop_label)
    plot_bandstructure(systems, phi, hoppings, title, 'v')

    return

def analyze_bandstructure_varyw():

    a      = 1
    b      = 1
    width  = 2
    phi    = 0.0
    u      = 1
    t_hops = [1.0]
    v_hops = [1.0]
    w_hops = [3.0, 2.0, 1.0, 0.0]
    
    hoppings = list((t,v,w) for t in t_hops for v in v_hops for w in w_hops)

    systems = []
    for hopping in hoppings:
        systems.append(make_system(a, b, phi, u, hopping, width).finalized())

    thop_label = round(hoppings[0][0], 3)
    vhop_label = round(hoppings[0][1], 3)
    title = 'Simple ladder band structure - t = {}eV, v = {}eV'.format(thop_label, vhop_label)
    plot_bandstructure(systems, phi, hoppings, title, 'w')

    return

def analyze_bandstructure_varyphi():

    a      = 1.0
    b      = 1.0
    width  = 2
    phi_steps = 8
    phis   = [2*pi*i/phi_steps for i in range(phi_steps)]
    u      = 5
    t_hops = 1
    v_hops = 1
    w_hops = 1

    hoppings = (t_hops, v_hops, w_hops)

    systems = []
    for phi in phis:
        systems.append(make_system(a, b, phi, u, hoppings, width).finalized())

    thop_label = round(hoppings[0], 3)
    vhop_label = round(hoppings[1], 3)
    whop_label = round(hoppings[2], 3)
    title = 'Simple ladder band structure - t = {}eV, v = {}eV, w = {}eV'.format(thop_label, vhop_label, whop_label)
    plot_bandstructure(systems, phis, hoppings, title, 'phi')

    return

def plot_spectrum(sys, n_eigs):

    colors = iter(cm.rainbow(np.linspace(0,1,n_eigs)))

    H_matrix = sys.hamiltonian_submatrix(sparse=True)
    eigenvals, eigenvecs = sorted_eigs(sla.eigsh(H_matrix.tocsc(), k=n_eigs, sigma=0))
    
    fig, ax = plt.subplots(figsize=(4,6))
    for n in range(n_eigs):
        ax.hlines(eigenvals[n], -1, 1, color=next(colors), label='n = ' + str(n))

    ax.legend()
    ax.set_title('Energy eigenvalues')
    ax.set_ylabel(r'$\epsilon$')
    ax.set_xlim(-1.5,3)
    ax.set_xticks([])

    plt.tight_layout()
    plt.show()

def plot_wavefunction(sys, n_eigs, plot_eigs):

    H_matrix = sys.hamiltonian_submatrix(sparse=True)
    eigenvals, eigenvecs = sorted_eigs(sla.eigsh(H_matrix.tocsc(), k=n_eigs, sigma=0))

    pos = range(1,len(plot_eigs)+1)

    fig = plt.figure(figsize=(8,12))
    for n in plot_eigs:
        ax = fig.add_subplot(len(plot_eigs), 1, pos[n])
        ax.plot(np.real(eigenvecs[::2,n]), label = r'$Re(\Psi)$')
        ax.plot(np.imag(eigenvecs[::2,n]), label = r'$Im(\Psi)$')

        ax.set_title(r'$\Psi_{state}(x,0)\quad \epsilon_{state}={energy}$'.format(state=n, energy=round(eigenvals[n],5)))
        ax.set_xlabel(r'$x\:[a_{}]$'.format(0))
        ax.set_ylabel(r'$\Psi$'.format(0))

    plt.tight_layout()
    plt.show()

    return

def plot_spec_prob(sys, n_eigs, plot_eigs):

    H_matrix = sys.hamiltonian_submatrix(sparse=True)
    eigenvals, eigenvecs = sorted_eigs(sla.eigsh(H_matrix.tocsc(), k=n_eigs, sigma=0))

    k = np.pi/4
    r = np.linspace(0, 100, 100)

    u_nk = []
    for i in range(n_eigs):
        u_nk.append(eigenvecs[::2, i]*np.exp(-1j*k*r))

    u_nk = np.array(u_nk)
    
    colors = iter(cm.rainbow(np.linspace(0,1,n_eigs)))
    pos = range(1,len(plot_eigs)+1)

    fig = plt.figure(figsize=(16,10))

    spec_ax = fig.add_subplot(1, 3, 1)
    for n in range(n_eigs):
        spec_ax.hlines(eigenvals[n], -1, 1, color=next(colors), label='n = ' + str(n))

    spec_ax.set_ylabel(r'$\epsilon\:[eV]$')
    spec_ax.set_xticks([])
    spec_ax.legend()

    for n in reversed(plot_eigs):
        wave_pos = 3*len(plot_eigs) - (pos[n] + 2*n)
        wave_ax = fig.add_subplot(len(plot_eigs), 3, wave_pos)
        wave_ax.plot(np.real(eigenvecs[::2, n]), label = r'$Re(\Psi)$')
        wave_ax.plot(np.imag(eigenvecs[::2, n]), label = r'$Im(\Psi)$')
        wave_ax.plot(np.imag(u_nk[n,:]), label = r'$Re(\Psi)$')
        wave_ax.set_title(r'$\Psi_{}(x,0)$'.format(n))
        wave_ax.set_xlabel(r'$x$')
        wave_ax.set_xlim(0,np.size(eigenvecs[::2,n])-1)

        dens_pos = 3*len(plot_eigs) - (pos[n]+ 2*n) + 1# Complicated due to the way positions are indexed in subplots. Could be an easier way?
        dens_ax = fig.add_subplot(len(plot_eigs), 3, dens_pos)
        kwant.plotter.map(sys, np.abs(eigenvecs[:,n])**2, colorbar=True, ax=dens_ax)
        dens_ax.set_title(r'$|\Psi_{state}(x,y)|^2\quad \epsilon_{state}={energy}$'.format(state=n, energy=round(eigenvals[n],5)))
        dens_ax.yaxis.tick_right()

    plt.tight_layout()
    plt.show()

    return       

def plot_band_state(sys):

    # What states to plot, momenta to choose from. Likely function parameters soon
    states = [0,1,2]
    state_momenta = [0,pi/4,pi/2,3*pi/4,pi]
    
    '''Retrieve positon/wavevector parameters'''
    # Extract site positions from system
    # Only the first half of the sites are in the 'fundamental domain'
    # Width is the number of y-values
    # Momenta in the first BZ are |k| < pi/b
    pos = np.array([site.pos for site in sys.sites])
    pos = pos[:len(pos)//2]     
    kx = state_momenta[4]


    '''Start the figure for band structure and eigenstates'''
    # Create band structure subplots
    # Add subplot in a 1 row, 2 column config. band_ax is axis 1 (counting from top left)
    # Initialize color interator (different color for each plotted eigenstate)
    fig = plt.figure(figsize=(10,int(len(states)/2)+6))
    band_ax = fig.add_subplot(1,2,1) 
    colors = iter(cm.rainbow(np.linspace(0,1,len(states))))


    '''Calculate band energies and momenta'''
    # Create a physic.Bands instance from the system
    # Define momenta for |k| < pi
    # Calculate eigenenergies by calling bands instance
    bands = kwant.physics.Bands(sys)
    momenta = np.linspace(0, pi, 100)
    energies = [bands(k) for k in momenta]


    '''Plot band structure on band_ax'''
    band_ax.plot(momenta, energies, color = 'k')
    band_ax.set_xlim(0,pi)
    band_ax.set_ylabel(r'$\epsilon(k)$')
    band_ax.set_xlabel(r'$k$')


    '''Plot range of eigenstates on state_ax'''
    # Calculate eigvals/vecs from bands instance 
    for n in states:

        # Retrieve eigvals/vecs from bands instance
        # eigvals = E_n(k_x)
        # eigvecs = u_{n,k_x}(0,y)
        eigvals, eigvecs  = bands(kx, return_eigenvectors=True)

        N = 10
        # Repeat u_nk in x-direction N times
        u_nk = []
        for _ in range(N): 
            u_nk.append(eigvecs[n,:])
        u_nk = np.array(u_nk)
        print('u_nk shape: ' + str(np.shape(u_nk)))

        # Calculate psi by modulating u_nk with plane wave
        psi = []
        for xi in range(N):
            psi.append(np.exp(1j*kx*xi)*u_nk[xi,:])
        psi = np.array(psi)
        print(r"$\Psi$ shape: " + str(np.shape(psi)))

        # Print some user output
        print("Momenta: kx = {:3.5f}".format(kx))

        # Iterate to next color for each plotted state
        color = next(colors)

        # Plot eigenstates on the state_ax
        # Mark the corresponding E_n(k) on band_ax
        # Plot \Psi_{nk}(y) on state_ax
        state_ax = fig.add_subplot(len(states), 2, 2*len(states) - 2*n)
        band_ax.scatter(kx, eigvals[n], color = color, zorder = 2, clip_on = False)
        state_ax.imshow(psi.real.transpose())
        state_ax.set_title('$\Psi_{}(x,y)$'.format(n), fontdict = {'color': color})
        state_ax.set_xlabel('x')
        state_ax.set_ylabel('y')
        state_ax.yaxis.tick_right()

    plt.show()

    return

def sorted_eigs(ev):
    evals, evecs = ev
    evals, evecs = map(np.array, zip(*sorted(zip(evals, evecs.transpose()))))
    return evals, evecs.transpose()

if __name__ == '__main__':
    main()
