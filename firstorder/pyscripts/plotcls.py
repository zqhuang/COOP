import numpy as np
from matplotlib import use
use('pdf')
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import pylab as m

fig, axes = plt.subplots(nrows=2, figsize=(8., 8.), sharex = True)

norm = (1./2.726e6/1.e-5)**2

cldata = np.loadtxt("recomb_slice_zetaCls.dat")
ells = [ s[0] for s in cldata ] 
clzz = [ s[3]*norm/5. for s in cldata ]
clzt = [ s[5]*norm for s in cldata ]
clze = [ s[6]*norm*5 for s in cldata ]

axes[0].plot(ells,clzz, color='r', label='$1/5 \zeta\zeta$')
axes[0].plot(ells,clzt, color='b', linestyle='--', label='$T\zeta$')
axes[0].plot(ells,clze, color='k', linestyle=':', label = '$5 E\zeta$')
axes[0].set_xscale('log')
axes[0].set_ylabel('$10^{10}\ell(\ell+1)C_\ell / (2\pi)$')
axes[0].set_xlabel('$\ell$')
axes[0].set_title('single slice')
axes[0].legend()

cldata = np.loadtxt("vis_zetaCls.dat")
ells = [ s[0] for s in cldata ] 
clzz = [ s[3]*norm/5. for s in cldata ] 
clzt = [ s[5]*norm for s in cldata ] 
clze = [ s[6]*norm*5 for s in cldata ] 

axes[1].plot(ells,clzz, color='r', label='$1/5\zeta\zeta$')
axes[1].plot(ells,clzt, color='b', linestyle='--', label='$T\zeta$')
axes[1].plot(ells,clze, color='k', linestyle=':', label = '$5 E\zeta$')
axes[1].set_xscale('log')
axes[1].set_ylabel('$10^{10}\ell(\ell+1)C_\ell / (2\pi)$')
axes[1].legend()
axes[1].set_title('visibility weighed')


plt.savefig('zetaCls.pdf', format='pdf')
