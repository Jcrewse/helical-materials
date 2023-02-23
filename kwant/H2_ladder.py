import kwant
from math import cos, exp

# Define the system Hamiltonian
# -----------------------------------------------------------------------------
def hamiltonian(a, b, width, hops, phi, N_phi = 1, finite = False):

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
    print(exp(1j*phi))
    l_v = (0.5*a**2)*(1-cos(phi)) + b**2 # Modulate hopping by geometric distance due to phi
    sys[kwant.builder.HoppingKind((1,0), lat, lat)]  = -v*exp(1j*phi)/l_v

    # Cross hoppings
    l_u = (0.5*a**2)*(1+cos(phi)) + b**2
    sys[kwant.builder.HoppingKind((1,1), lat, lat)]  = -w*exp(1j*phi)/l_u
    sys[kwant.builder.HoppingKind((-1,1), lat, lat)] = -w*exp(1j*phi)/l_u

    return sys