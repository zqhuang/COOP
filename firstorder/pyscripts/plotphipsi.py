import numpy as np
from matplotlib import use
use('pdf')
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.ylabel('$\Phi/\Psi$')
plt.xlabel('$a$')
plt.xscale('log')

data0 = np.loadtxt("pert_alphaM0.txt")
datap = np.loadtxt("pert_alphaM1.txt")
datan = np.loadtxt("pert_alphaM-1.txt")

a0 = np.array([ np.exp(s[0]) for s in data0 ])
ap = np.array([ np.exp(s[0]) for s in datap ])
an = np.array([ np.exp(s[0]) for s in datan ])
rat0 = np.array([ s[5]/s[6]-1. for s in data0 ])
ratp = np.array([ s[5]/s[6]-1. for s in datap ])
ratn = np.array([ s[5]/s[6]-1. for s in datan ])

plt.plot(a0, rat0, color='k')
plt.plot(ap, ratp, color='g')
plt.plot(an, ratn, color='b')
plt.savefig('PhibyPsi_alphaM.pdf', format='pdf')
