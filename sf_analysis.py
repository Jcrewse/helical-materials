import pickle
import numpy as np
import itertools
import matplotlib.pyplot as plt
from gpaw.unfold import make_colormap

emin = -15
emax = 15

npts = 10000
width = 0.0002

colors = ['blue', 'red', 'green']

Nk = 100#int(input('N k-points: '))

plt.figure()
colors_tuple = ('blue', 'red', 'green')
colors = itertools.cycle(colors_tuple)

Nphis = [4,5,8]

for i in Nphis:
    
    system_name = f'H2-Nphi-{i}'#input('System tag: ')

    filename = f'sf_{system_name}_unfolded'
    e, A_ke, x, X, points_name = pickle.load(open(filename + '.pckl','rb'))
    print('--------{}'.format(filename))
    print(x)

    A_ke /= np.max(A_ke)
    A_ek = A_ke.T
    A_ekc = np.reshape(A_ek, (A_ek.shape[0], A_ek.shape[1]))

    color = next(colors)
    mycmap = make_colormap(color)

    #plt.plot([0, x[-1]], 2 * [0.0], '--', c='0.5')
    plt.imshow(A_ekc+0.23,
                cmap=mycmap,
                aspect='auto',
                origin='lower',
                vmin=0.,
                vmax=1,
                extent=[0, Nphis[-1]*np.pi, e.min(), e.max()])

    # for k in X[1:-1]:
    #     plt.plot([k, k], [emin, emax], lw=0.5, c='0.5')

plt.xticks([0,Nphis[-1]*np.pi], points_name, size=20)
plt.ylabel('E(eV)')
plt.axis([0, Nphis[-1]*np.pi, emin, emax])
plt.tight_layout()
plt.savefig('sf.png')
plt.show()



