import numpy as np
from matplotlib import use
use('pdf')
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import pylab as m

cdict = {
'red'  :  ((0., 1., 1.), (0.03, 0.56, 0.56), (0.2, 0., 0.), (0.3, 0., 0.), (0.4, 0.1, 0.1), (0.5, 0., 0.) , (0.6, 0.7, 0.7), (0.7, 1., 1.), (0.8, 1., 1.), (1., 1., 1.)),
'green':  ((0., 1., 1.), (0.03, 0.21, 0.21), (0.2, 0., 0.), (0.3, 1., 1.), (0.4, 0.9, 0.9) ,(0.5, 1., 1.),  (0.7, 1., 1.),  (0.8, 0.5, 0.5), (1., 0., 0.)),
'blue' :  ((0., 1., 1.), (0.03, 0.94, 0.94), (0.2, 1., 1.), (0.3, 1., 1.), (0.4, 0.1, 0.1), (0.5, 0., 0.),  (0.7, 0., 0.), (0.8, 0.2, 0.2),  (1., 0., 0.))
}
#generate the colormap with 1024 interpolated values
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
mpl.cm.register_cmap(cmap = my_cmap)

plt.set_cmap(my_cmap) #'gist_rainbow_r')

fig, axes = plt.subplots(nrows=4, figsize=(8., 8.), sharex = True)


chimin = 0.10000E-01
chimax = 1.2000
lmin = 2.
lmax = 3000.
minsn = 0.01
maxsn = 10.
chi_z10 = 0.694
chi_z1 = 0.245
chi_zhalf = 0.141

axes[0].set_yscale('log')
axes[0].set_ylabel('$\ell$')
grid = np.loadtxt("nztinfo.txt").transpose()
im = axes[0].imshow(grid[::-1], extent=[chimin,chimax,lmin, lmax], vmin = minsn, vmax = maxsn, norm=LogNorm())
axes[0].set_title('Temperature $S/N$')

axes[1].set_yscale('log')
axes[1].set_ylabel('$\ell$')
grid = np.loadtxt("nzeinfo.txt").transpose()
im = axes[1].imshow(grid[::-1], extent=[chimin,chimax,lmin, lmax], vmin = minsn, vmax = maxsn, norm=LogNorm())
axes[1].set_title('E polarization $S/N$')

axes[2].set_yscale('log')
axes[2].set_ylabel('$\ell$')
grid = np.loadtxt("nzetinfo.txt").transpose()
im = axes[2].imshow(grid[::-1], extent=[chimin,chimax,lmin, lmax], vmin = minsn, vmax = maxsn, norm=LogNorm())
axes[2].set_title('T + E $S/N$')


fig.colorbar(im, ax = axes.ravel().tolist())

chivis = np.loadtxt("vis.txt")
chi = [ s[0] for s in chivis ]
vis = [ s[1] for s in chivis ]
axes[3].plot(chi,vis)
axes[3].set_yscale('log')
axes[3].text(0.1, 1., "ISW")
axes[3].axvline(x = chi_z1, color="black", ls="--")
axes[3].text(chi_z1 + 0.02, 1., '$z=1$', rotation = 90)
axes[3].text(0.58, 1., "reionization")
axes[3].text(0.9, 1., "recombination")


axes[3].set_ylabel('$\dot\kappa e^{-\kappa} / H_0$')
axes[3].set_xlabel('$\chi/\chi_{rec}$')
axes[3].set_title('differential visibility')

plt.savefig('zetaSN_with_Noise.pdf', format='pdf')
