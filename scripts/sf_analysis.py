import pickle
import numpy as np
import itertools
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from gpaw.unfold import make_colormap

emin = -35
emax = 20

npts = 10000
width = 0.0002

Nk = 100#int(input('N k-points: '))

plt.figure(figsize=(6,8))
colors_tuple = ('Blues', 'Reds', 'Greens', 'Purples', 'Oranges')
colors = itertools.cycle(colors_tuple)

Nphis = [1,2,3,4,6]

for i in Nphis:
    
    system_name = f'CO-Nphi-{i}'#input('System tag: ')

    # 
    filename = f'sf_{system_name}_unfolded'
    e, A_ke, x, X, points_name = pickle.load(open(filename + '.pckl','rb'))
    print('--------{}'.format(filename))

    # Normalize spectral functions and transpose for plotting
    A_ke /= np.max(A_ke)
    A_ek = A_ke.T
    A_ekc = np.reshape(A_ek, (A_ek.shape[0], A_ek.shape[1]))

    # Cycle to the next color scheme
    color = next(colors)
    alphas = np.zeros_like(A_ekc)
    alphas[A_ekc > 0.01] = 1

    #plt.plot([0, x[-1]], 2 * [0.0], '--', c='0.5')
    plt.imshow(A_ekc+0.5,
                cmap=color,
                aspect='auto',
                origin='lower',
                vmin=0.,
                vmax=1,
                alpha=alphas,
                label=r'$N_{\phi} = $' + f'{i}',
                extent=[0, Nphis[-1]*np.pi, e.min(), e.max()])

plt.xticks([0,Nphis[-1]*np.pi], points_name, size=20)
plt.ylabel('E(eV)')
plt.axis([0, Nphis[-1]*np.pi, emin, emax])
plt.legend()
plt.tight_layout()
plt.savefig('sf.png')
plt.show()



