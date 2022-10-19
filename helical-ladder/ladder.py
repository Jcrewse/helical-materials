import kwant
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import scipy.sparse.linalg as sla
from math import pi, cos

# Conversion factors FROM SI TO ATOMIC
# Parameters in the code are written in SI (energy = eV, length=angstrom, etc...)
# Multiply written parameter by these to get atomic units (hbar = m = c = 1)
CONVERT_HARTREE = (27.211386)**-1
CONVERT_EV      = CONVERT_HARTREE**-1
CONVERT_BOHR_R  = (0.529177)**-1
CONVERT_ANGS    = CONVERT_BOHR_R**-1

hop_dict = {'t':0, 'v':1, 'w':2}

def main():

    analyze_bandstructure_varyt()
    analyze_bandstructure_varyv()
    analyze_bandstructure_varyw()
    analyze_bandstructure_varyphi()

    # Create a finite system to visualize the structure
    system = make_system(a=1, b=1, phi=np.pi/4, u=5, hops=(1,1,1), width=2, length=50).finalized()
    #plot_structure(system)
    #plot_spectrum(system, n_eigs=10)
    #plot_wavefunction(system, n_eigs=10, plot_eigs=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    plot_spec_wave(system, n_eigs=10, plot_eigs=[0,1,2,3,4,5])
    #translation(system)

    return

def make_system(a, b, phi, u, hops, width, length=None):

    t, v, w = hops

    # Perform the calculation in atomic units
    u = u*CONVERT_HARTREE
    t = t*CONVERT_HARTREE
    v = v*CONVERT_HARTREE
    w = w*CONVERT_HARTREE

    a = a*CONVERT_BOHR_R
    b = b*CONVERT_BOHR_R


    lat = kwant.lattice.general([(b,0), (0,a)], norbs=1)

    sys = kwant.Builder()

    if length is not None:
        # Create a finite sizes system to show structure
        sys[lat.shape((lambda pos: (pos[1] >= 0 and pos[1] < a*width) and (pos[0] >= 0 and pos[0] < b*length)), (0,0))] = u
    else:
        # Create an infinite system for band structure calculations
        sys = kwant.Builder(kwant.TranslationalSymmetry((b,0)))
        sys[lat.shape((lambda pos: pos[1] >= 0 and pos[1] < a*width), (0,0))] = u

    # Intra-layer hoppings
    sys[kwant.builder.HoppingKind((0,1), lat, lat)]  = -t

    # Inter-layer hoppings
    sys[kwant.builder.HoppingKind((1,0), lat, lat)]  = -v*((0.5*a**2)*(1-cos(phi)) + b**2)**(-1)

    # Cross hoppings
    sys[kwant.builder.HoppingKind((1,1), lat, lat)]  = -w*((0.5*a**2)*(1+cos(phi)) + b**2)**(-1)
    sys[kwant.builder.HoppingKind((-1,1), lat, lat)] = -w*((0.5*a**2)*(1+cos(phi)) + b**2)**(-1)

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
        kwant.plotter.map(sys, np.abs(eigenvecs[:,n])**2, colorbar=True, ax=ax)

        ax.set_title(r'$|\Psi_{state}(x,y)|^2\quad \epsilon_{state}={energy}$'.format(state=n, energy=round(eigenvals[n],5)))
        ax.set_xlabel(r'$x\:[a_{}]$'.format(0))
        ax.set_ylabel(r'$y\:[a_{}]$'.format(0))

    plt.show()

    return

def plot_spec_wave(sys, n_eigs, plot_eigs):

    H_matrix = sys.hamiltonian_submatrix(sparse=True)
    eigenvals, eigenvecs = sorted_eigs(sla.eigsh(H_matrix.tocsc(), k=n_eigs, sigma=0))
    
    colors = iter(cm.rainbow(np.linspace(0,1,n_eigs)))
    pos = range(1,len(plot_eigs)+1)

    fig = plt.figure(figsize=(16,8))

    spec_ax = fig.add_subplot(1, 2, 1)
    for n in range(n_eigs):
        spec_ax.hlines(eigenvals[n], -1, 1, color=next(colors), label='n = ' + str(n))

    spec_ax.set_ylabel(r'$\epsilon\:[eV]$')
    spec_ax.set_xticks([])
    spec_ax.legend()

    for n in reversed(plot_eigs):
        wave_pos = 2*len(plot_eigs) - (pos[n] + n) + 1 # Complicated due to the way positions are indexed in subplots. Could be an easier way?
        wave_ax = fig.add_subplot(len(plot_eigs), 2, wave_pos)
        kwant.plotter.map(sys, np.abs(eigenvecs[:,n])**2, colorbar=True, ax=wave_ax)

        wave_ax.set_title(r'$|\Psi_{state}(x,y)|^2\quad \epsilon_{state}={energy}$'.format(state=n, energy=round(eigenvals[n],5)))
        wave_ax.yaxis.tick_right()

    plt.show()

    return

def translation(sys):

    H_matrix = sys.hamiltonian_submatrix(sparse=True)
    eigenvals, eigenvecs = sorted_eigs(sla.eigsh(H_matrix.tocsc(), k=5, sigma=0))
    print(np.shape(eigenvecs))

    return            

def sorted_eigs(ev):
    evals, evecs = ev
    evals, evecs = map(np.array, zip(*sorted(zip(evals, evecs.transpose()))))
    return evals, evecs.transpose()

if __name__ == '__main__':
    main()
