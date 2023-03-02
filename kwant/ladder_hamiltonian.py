import kwant
import numpy as np
from math import cos

# Define the system Hamiltonian
# -----------------------------------------------------------------------------
def hamiltonian(a, b, width, hops, phi, N_phi = 1, finite = False):

    # Unpack the hoppings tuple
    u, t, v, w = hops

    # Define the lattice
    lat_vecs = [(a,0), (0,b)]
    lat = kwant.lattice.general(lat_vecs, norbs=1)

    # Define Hamiltonian
    sys = kwant.Builder()
    if finite:
        # Create a finite sizes system to show structure
        sys[lat.shape((lambda pos: (pos[1] >= 0 and pos[1] < b*width) 
                       and (pos[0] >= 0 and pos[0] < N_phi*a)), (0,0))] = u
    else:
        # Create an infinite system for band structure calculations
        sys = kwant.Builder(kwant.TranslationalSymmetry((N_phi*a,0)))
        sys[lat.shape((lambda pos: (pos[1] >= 0 and pos[1] < b*width) 
                       and (pos[0] >= 0 and pos[0] < N_phi*a)), (0,0))] = u

    # Intra-layer hoppings
    sys[kwant.builder.HoppingKind((0,1), lat, lat)]  = -t

    # Inter-layer hoppings
    l_v = (0.5*a**2)*(1-np.cos(phi)) + b**2 # Modulate hopping by geometric distance due to phi
    sys[kwant.builder.HoppingKind((1,0), lat, lat)]  = -v*np.exp(-1j*phi)/l_v

    # Cross hoppings
    l_u = (0.5*a**2)*(1+np.cos(phi)) + b**2
    sys[kwant.builder.HoppingKind((1,1), lat, lat)]  = -w*np.exp(-1j*phi)/l_u
    sys[kwant.builder.HoppingKind((-1,1), lat, lat)] = -w*np.exp(-1j*phi)/l_u

    # Output 
    print('========== Hamiltonian initialized ==========')
    print('Lattice constants: a = {}, b = {}'.format(a, b))
    print('Finite width = {}'.format(width))
    print('Hoppings: u = {}, t = {}, v = {}, w = {}'.format(u,t,v,w))
    if (N_phi==1): print('Twist angle: phi = 0')
    if (N_phi>1): print('Twist angle: phi = {}'.format(round(phi,3)))
    print('---------- Hamiltonian matrix ----------')
    print('H = \n' + str(sys.finalized().hamiltonian_submatrix()))
    print('H_cell = \n' + str(sys.finalized().cell_hamiltonian()))
    print('H_inter_cell = \n' + str(sys.finalized().inter_cell_hopping()))

    return sys