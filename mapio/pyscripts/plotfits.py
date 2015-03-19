import numpy as np
from matplotlib import use
use('pdf')
import matplotlib.pyplot as plt
import healpy as hp


m = hp.read_map("../planck14/smica_masked.fits")
m = m * 1.e6
hp.mollview(map = m, title = "Planck Full-mission Temperature Map", cbar = True, min = -420., max = 420., unit="$\mu K$", xsize = 600)
plt.savefig("planckmap_full_mission.pdf", format="pdf")
