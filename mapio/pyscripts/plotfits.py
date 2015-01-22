import numpy as np
from matplotlib import use
use('pdf')
import matplotlib.pyplot as plt
import healpy as hp


NSIDE = 512
m = hp.read_map("../dust_peaks.fits")
hp.mollview(map = m, title = "peak catalog $e>0.6,\ I>0.5 \sigma_0$", cbar = False)
plt.savefig("peaks.pdf", format="pdf")
