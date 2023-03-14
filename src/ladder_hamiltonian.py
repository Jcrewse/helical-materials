import numpy as np
import kwant
from math import cos, sin

def hamiltonian(lat_const_a = 1, lat_const_b = 1, width = 2, hops = (20,25,10,5), phi = 0,
                n_phi = 1, output=False):
    """Returns the kwant.Builder object of the helical ladder system with the
    twist angle implemented as a unitary transformation of the molecular 
    Hamtiltonian."""

    # Unpack the hoppings tuple
    on_site_pot, intra_hop, inter_hop, cross_hop = hops

    # Define the lattice
    lat_vecs = [(lat_const_a,0), (0,lat_const_b)]
    lat = kwant.lattice.general(lat_vecs, norbs=1)

    # System is translationally symmetric for translations
    # T = n_phi*a
    symm = kwant.lattice.TranslationalSymmetry((n_phi*lat_const_a,0))
    sys  = kwant.builder.Builder(symm)
    
    # Define the molecular Hamiltonian
    H_molecule = [[on_site_pot, -intra_hop],[-intra_hop, on_site_pot]]
    
    # Rotation matrices
    rot_matrix = [[cos(phi), -sin(phi)],[sin(phi), cos(phi)]]
    rot_matrix_dagger = np.transpose(np.conjugate(rot_matrix))
    
    H_rotated = np.array(H_molecule)
    for n in range(n_phi):
        # Rotate the molecular Hamiltonian
        # For each layer in unit cell, perform the unitary transformation
        # H --> RHR^+
        # This transformation is applied each layer so that after 
        # n_phi layers we have
        # H --> R^{n_phi}HR^{+ n_phi}
        
        # Carry out the unitary transformation once per layer
        if n >= 1:
            H_rotated = np.matmul(H_rotated, rot_matrix_dagger)
            H_rotated = np.matmul(rot_matrix, H_rotated)
        
        # On-site pot in each layer set to diagonal terms 
        sys[lat(n,0)] = H_rotated[0,0]
        sys[lat(n,1)] = H_rotated[1,1]
        
        # Intra-layer hoppings in each layer set as off-diagonals
        sys[lat(n,0), lat(n,1)] = H_rotated[0,1]
        sys[lat(n,1), lat(n,0)] = H_rotated[1,0]
    
    # Inter-layer hoppings
    l_v = 1#(0.5*lat_const_a**2)*(1-np.cos(phi)) + lat_const_b**2
    sys[kwant.builder.HoppingKind((1,0), lat, lat)]   = -inter_hop/l_v

    # Cross hoppings
    l_u = 1#(0.5*lat_const_a**2)*(1+np.cos(phi)) + lat_const_b**2
    sys[kwant.builder.HoppingKind((1,1), lat, lat)]   = -cross_hop/l_u
    sys[kwant.builder.HoppingKind((-1,1), lat, lat)]  = -cross_hop/l_u
    
    # Output
    if output:
        np.set_printoptions(linewidth=500)
        print('========== Hamiltonian initialized ==========')
        print(f'Lattice constants: a = {lat_const_a}, b = {lat_const_b}')
        print(f'Finite width = {width}')
        print(f'Hoppings: u = {on_site_pot}, t = {intra_hop}, v = {inter_hop}, w = {cross_hop}')
        if n_phi==1:
            print('Twist angle: phi = 0')
        if n_phi>1:
            print(f'Twist angle: phi = {phi}')
        print('---------- Hamiltonian matrix ----------')
        print('H = \n' + str(sys.finalized().hamiltonian_submatrix()))
        print('H_cell = \n' + str(sys.finalized().cell_hamiltonian()))
        print('H_inter_cell = \n' + str(sys.finalized().inter_cell_hopping()))

    return sys

###############################################################################
###############################################################################
###############################################################################

def old_hamiltonian(lat_const_a, lat_const_b, width, hops, phi,
                n_phi = 1, finite = False, output=False):
    """Returns the kwant.Builder object of the helical ladder system based with
    twist angle implemented as simple geometric modulation of the hoppings."""

    # Unpack the hoppings tuple
    on_site_pot, intra_hop, inter_hop, cross_hop = hops

    # Define the lattice
    lat_vecs = [(lat_const_a,0), (0,lat_const_b)]
    lat = kwant.lattice.general(lat_vecs, norbs=1)

    # Define Hamiltonian
    sys = kwant.Builder
    if finite:
        sys[lat.shape((lambda pos: (pos[1] >= 0 and pos[1] < lat_const_b*width) 
                       and (pos[0] >= 0 and pos[0] < n_phi*lat_const_a)), (0,0))] = on_site_pot
    else:
        # System is translationally symmetric for n_phi layers
        sys = kwant.builder.Builder(kwant.lattice.TranslationalSymmetry((n_phi*lat_const_a,0)))
        
        # On-site potential the same for each lattie site on the ladder
        sys[lat.shape((lambda pos: (pos[1] >= 0 and pos[1] < lat_const_b*width)
                       and (pos[0] >= 0 and pos[0] < n_phi*lat_const_a)), (0,0))] = on_site_pot

    # Intra-layer hoppings
    sys[kwant.builder.HoppingKind((0,1), lat, lat)]  = -intra_hop

    # Inter-layer hoppings
    l_v = (0.5*lat_const_a**2)*(1-np.cos(phi)) + lat_const_b**2
    sys[kwant.builder.HoppingKind((1,0), lat, lat)]   = -inter_hop/l_v

    # Cross hoppings
    l_u = (0.5*lat_const_a**2)*(1+np.cos(phi)) + lat_const_b**2
    sys[kwant.builder.HoppingKind((1,1), lat, lat)]   = -cross_hop/l_u
    sys[kwant.builder.HoppingKind((-1,1), lat, lat)]  = -cross_hop/l_u
    
    # Output
    if output:
        print('========== Hamiltonian initialized ==========')
        print(f'Lattice constants: a = {lat_const_a}, b = {lat_const_b}')
        print(f'Finite width = {width}')
        print(f'Hoppings: u = {on_site_pot}, t = {intra_hop}, v = {inter_hop}, w = {cross_hop}')
        if n_phi==1:
            print('Twist angle: phi = 0')
        if n_phi>1:
            print(f'Twist angle: phi = {phi}')
        print('---------- Hamiltonian matrix ----------')
        print('H = \n' + str(np.matrix(sys.finalized().hamiltonian_submatrix())))
        print('H_cell = \n' + str(sys.finalized().cell_hamiltonian()))
        print('H_inter_cell = \n' + str(sys.finalized().inter_cell_hopping()))

    return sys
