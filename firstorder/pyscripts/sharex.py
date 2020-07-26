import numpy as np
from matplotlib import use
use('pdf')
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import pylab as m

cdict = {
'red'  :  ((0., 1., 1.), (0.02, 0.56, 0.56), (0.25, 0., 0.), (0.45, 0.1, 0.1), (0.5, 0., 0.) , (0.6, 0.7, 0.7), (0.7, 1., 1.), (0.8, 1., 1.), (1., 1., 1.)),
'green':  ((0., 1., 1.), (0.02, 0.21, 0.21), (0.25, 0., 0.), (0.45, 0.9, 0.9) ,(0.5, 1., 1.),  (0.7, 1., 1.),  (0.8, 0.5, 0.5), (1., 0., 0.)),
'blue' :  ((0., 1., 1.), (0.02, 0.94, 0.94), (0.25, 1., 1.), (0.45, 0.1, 0.1), (0.5, 0., 0.),  (0.7, 0., 0.), (0.8, 0.2, 0.2),  (1., 0., 0.))
}
#generate the colormap with 1024 interpolated values
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
mpl.cm.register_cmap(cmap = my_cmap)

plt.set_cmap(my_cmap) #'gist_rainbow_r')

fig, axes = plt.subplots(nrows=4, figsize=(10, 8.), sharex = True)


chimin = 0.10000E-01
chimax = 1.2000
lmin = 2.
lmax = 3000.

axes[0].set_yscale('log')
axes[0].set_ylabel('$\ell$')
grid = np.loadtxt("ztinfo.txt").transpose()
im = axes[0].imshow(grid[::-1], extent=[chimin,chimax,lmin, lmax], vmin = 0.1, vmax = 10., norm=LogNorm())
axes[0].set_title('Temperature $S/N$')

axes[1].set_yscale('log')
axes[1].set_ylabel('$\ell$')
grid = np.loadtxt("zeinfo.txt").transpose()
im = axes[1].imshow(grid[::-1], extent=[chimin,chimax,lmin, lmax], vmin = 0.1, vmax = 10., norm=LogNorm())
axes[1].set_title('E polarization $S/N$')

axes[2].set_yscale('log')
axes[2].set_ylabel('$\ell$')
grid = np.loadtxt("zetinfo.txt").transpose()
im = axes[2].imshow(grid[::-1], extent=[chimin,chimax,lmin, lmax], vmin = 0.1, vmax = 10., norm=LogNorm())
axes[2].set_title('T + E $S/N$')


fig.colorbar(im, ax = axes.ravel().tolist())

chivis = np.loadtxt("vis.txt")
chi = [ s[0] for s in chivis ]
vis = [ s[1] for s in chivis ]
axes[3].plot(chi,vis)
axes[3].set_yscale('log')
axes[3].set_ylabel('$\dot\kappa e^{-\kappa} / H_0$')
axes[3].set_xlabel('$\chi/\chi_{rec}$')
axes[3].set_title('differential visiblity')

plt.savefig('zetaSN.pdf', format='pdf')
