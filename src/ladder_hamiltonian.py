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
        # System is translationally symmetric for n_phi layers
        sys = kwant.builder.Builder(kwant.lattice.TranslationalSymmetry((n_phi*lat_const_a,0)))
        
        # H_11
        H_11 = cos(phi)*(intra_hop*sin(phi) + on_site_pot*cos(phi)) \
                        + sin(phi)*(intra_hop*cos(phi) + on_site_pot*sin(phi))
        H_22 = cos(phi)*(on_site_pot*cos(phi) - intra_hop*sin(phi)) \
                        + sin(phi)*(on_site_pot*sin(phi) - intra_hop*cos(phi))
        sys[lat(0,0)] = H_11
        sys[lat(0,1)] = H_22

    # Intra-layer hoppings
    # H_12
    H_12 = sin(phi)*(intra_hop*sin(phi) + on_site_pot*cos(phi)) \
            - cos(phi)*(intra_hop*cos(phi) + on_site_pot*sin(phi))
    sys[kwant.builder.HoppingKind((0,1), lat, lat)] = H_12
            
    # H_21
    H_21 = cos(phi)*(on_site_pot*sin(phi) - intra_hop*cos(phi)) \
            - sin(phi)*(on_site_pot*cos(phi) - intra_hop*sin(phi))
    sys[kwant.builder.HoppingKind((0,-1), lat, lat)] = H_21
    
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
