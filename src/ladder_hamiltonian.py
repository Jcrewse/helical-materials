import numpy as np
import kwant

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
        sys = kwant.builder.Builder(kwant.lattice.TranslationalSymmetry((n_phi*lat_const_a,0)))
        sys[lat.shape((lambda pos: (pos[1] >= 0 and pos[1] < lat_const_b*width)
                       and (pos[0] >= 0 and pos[0] < n_phi*lat_const_a)), (0,0))] = on_site_pot

    # Intra-layer hoppings
    sys[kwant.builder.HoppingKind((0,1), lat, lat)]  = -intra_hop

    # Inter-layer hoppings
    l_v = (0.5*lat_const_a**2)*(1-np.cos(phi)) + lat_const_b**2
    sys[kwant.builder.HoppingKind((1,0), lat, lat)]   = -inter_hop*np.exp(-1j*phi)/l_v

    # Cross hoppings
    l_u = (0.5*lat_const_a**2)*(1+np.cos(phi)) + lat_const_b**2
    sys[kwant.builder.HoppingKind((1,1), lat, lat)]   = -cross_hop*np.exp(-1j*(phi))/l_u
    sys[kwant.builder.HoppingKind((-1,1), lat, lat)]  = -cross_hop*np.exp(1j*(phi))/l_u
    
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
