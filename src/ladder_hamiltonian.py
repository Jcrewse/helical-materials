import numpy as np
import kwant
from math import cos, sin

# Define the system Hamiltonian
# -----------------------------------------------------------------------------
def hamiltonian(lat_const_a, lat_const_b, width, hops, phi,
                n_phi = 1, finite = False, output=False):
    """Returns the kwant.Builder object of the helical ladder system."""

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
        # System is translationally symmetric for translations
        # T = n_phi*a
        sys = kwant.builder.Builder(kwant.lattice.TranslationalSymmetry((n_phi*lat_const_a,0)))
        
        # Define the molecular Hamiltonian
        # Rotation matrices
        H_molecule = [[on_site_pot, -intra_hop],[-intra_hop, on_site_pot]]
        rot_matrix = [[cos(phi), -sin(phi)],[sin(phi), cos(phi)]]
        rot_matrix_dagger = np.transpose(rot_matrix)
        
        H_rotated = H_molecule
        for n in range(n_phi):
            # Rotate the molecular Hamiltonian
            # For each layer in unit cell, perform the unitary transformation
            # H --> RHR^+
            # This transformation is applied each layer so that after 
            # n_phi layers we have
            # H --> R^{n_phi}HR^{+ n_phi}
            print(f'Layer = {n}')
            print(f'H_molecule = {H_molecule}')
            # Raise rotation matrices to appropriate power
            R  = np.linalg.matrix_power(rot_matrix, n)
            RT = np.linalg.matrix_power(rot_matrix_dagger, n) 
            
            # Carry out the unitary transformation
            H_rotated = np.matmul(H_rotated, RT)
            H_rotated = np.matmul(R, H_rotated)
            print(f'H_rot = \n {H_rotated}')
            
            # On-site pot in each layer set to diagonal terms 
            sys[lat(n,0)] = H_rotated[0,0]
            sys[lat(n,1)] = H_rotated[1,1]

            # Intra-layer hoppings in each layer set as off-diagonals
            sys[kwant.builder.HoppingKind((n,1), lat, lat)] = H_rotated[0,1]
            sys[kwant.builder.HoppingKind((n,-1), lat, lat)] = H_rotated[1,0]
    
    # Inter-layer hoppings
    l_v = 1#(0.5*lat_const_a**2)*(1-np.cos(phi)) + lat_const_b**2
    sys[kwant.builder.HoppingKind((1,0), lat, lat)]   = -inter_hop/l_v

    # Cross hoppings
    l_u = 1#(0.5*lat_const_a**2)*(1+np.cos(phi)) + lat_const_b**2
    sys[kwant.builder.HoppingKind((1,1), lat, lat)]   = -cross_hop/l_u
    sys[kwant.builder.HoppingKind((-1,1), lat, lat)]  = -cross_hop/l_u
    
    # Output
    if output:
        np.set_printoptions(linewidth=200)
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
