import numpy as np
from matplotlib import use
use('pdf')
import matplotlib.pyplot as plt
import healpy as hp


## the path to the maps
mapdir = "./"
#"../planck14/"

##load the temperature map
tmap = hp.read_map(mapdir + "dx11_v2_smica_int_cmb_020a_0512.fits")
## unit conversion K => mu K 
tmap  = tmap *  1.e6   

##load mask
maskmap = hp.read_map(mapdir + "dx11_v2_common_int_mask_020a_0512.fits")

##apply the mask to the map
tmap_masked = hp.ma(tmap)
tmap_masked.mask = np.logical_not(maskmap)

##plot the masked map
fig1 = plt.figure()
hp.mollview(map = tmap_masked, title = "Planck Full-mission Temperature Map", cbar = True, min = -420., max = 420., unit="$\mu K$", xsize = 600)
plt.savefig("plancktmap.pdf", format="pdf")

##compute C_l's (using the inpainted fullsky map)
lmax = 20
twopi = 2.*np.pi
##compute beam width
arcmin = 1./60.*np.pi/180.
fwhm = 20.*arcmin
sigma = fwhm/np.sqrt(8.*np.log(2.))
def beam_function(l):
    """ the beam window function W(l) that smoothes the map (C_l is multiplied by W(l)^2)"""
    return np.exp(-l*(l+1.)*sigma**2/2.)

##compute Cl's
cls = hp.sphtfunc.anafast(tmap, lmax = lmax, alm = True)

#ells = []
#normed_cls = []
#for l in range(2, lmax+1):
#    ells.append( l*1.)
#    normed_cls.append ( l*(l+1.)*cls[l]/twopi / beam_function(l)**2 )
##plot Cl's    
#fig2 = plt.figure()    
#plt.plot(ells, normed_cls)
#plt.xlabel("$\ell$")
#plt.ylabel("$\ell(\ell+1)C_\ell/(2\pi) [\mu K^2]$")
#plt.savefig("Cls.pdf", format = "pdf")
    
