import numpy as np
from matplotlib import use
use('pdf')
import matplotlib.pyplot as plt
import healpy as hp
import sys
map = hp.read_map(sys.argv[1])
fig = plt.figure()
hp.mollview(map = map, title = sys.argv[3], cbar = True, xsize = 600)
plt.savefig(sys.argv[2], format="pdf")

