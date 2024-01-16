#!/usr/bin/env python

#!/usr/bin/env python

import numpy as np
import healpy as hp
from newsetup_matplotlib import *
from planckcolors import planck_parchment_cmap, planck_grey_cmap,colombi1_cmap
from matplotlib import cm
from plot import *
import idlsave
import pyfits as py

nside=512
width=18.0
cmap = planck_parchment_cmap
cmap2 = planck_grey_cmap
xsize = 2000
ysize = xsize/2.0


theta = np.linspace(np.pi, 0, ysize)
phi = np.linspace(-np.pi, np.pi, xsize)
longitude = np.radians(np.linspace(-180, 180, xsize))
latitude = np.radians(np.linspace(-90, 90, ysize))

# project the map to a rectangular matrix xsize x ysize
PHI, THETA = np.meshgrid(phi, theta)
grid_pix = hp.ang2pix(nside, THETA, PHI)

from matplotlib.projections.geo import GeoAxes

class ThetaFormatterShiftPi(GeoAxes.ThetaFormatter):
    """Shifts labelling by pi

    Shifts labelling from -180,180 to 0-360"""
    def __call__(self, x, pos=None):
        if x != 0:
            x *= -1
        if x < 0:
            x += 2*np.pi
        return GeoAxes.ThetaFormatter.__call__(self, x, pos)

nrow=3
ncol=2
rat = float(nrow)/float(ncol)

height = 0.6 * width * rat
fig = plt.figure(figsize=(cm2inch(width), cm2inch(height)))

count=0
title=[r'$D_{353}$',r'$D_{353}^{\mathrm {b}}$', r'$\lambda_{-}$', r'$\lambda_{+}$', r'filaments $\lambda_{-}$', r'filaments $\lambda_{+}$']

fits1 = py.open('fits/Planck_dustmodel_353GHz_15arcmin_ns512.fits')
model=fits1[0].data

fil_model = hp.read_map('fits/input_dustmap_n2_nscale5.fits',0)

m1=idlsave.read('savesets/lambda_minus_plus.sav')

s=idlsave.read('savesets/filaments_data_n2_nside_512_rotangle_15_cutpixel_20_datatype_filter.sav')
map1=m1.lambda_minus
map2=s.filament_map
mm7=s.length_map

map2[np.where(mm7 < 2.1)]=0.0

s2=idlsave.read('savesets/plus_filaments_data_n2_nside_512_rotangle_15_cutpixel_20_datatype_filter.sav')
map3=m1.lambda_plus
map4=s2.filament_map_plus
mm8=s2.length_map_plus

map4[np.where(mm8 <= 2.1)]=0.0

mask= hp.read_map("fits/mask_gal_cut_30.0000_LMC_SMC.fits",0)

for i in range(6):


    if(i == 0): m = model[:]*1.e6
    if(i == 1): m = fil_model[:]*1.e6
    if(i == 2): m = map1[:]
    if(i == 3): m = map3[:]
    if(i == 4): m = map2[:]
    if(i == 5): m = map4[:]
 
 #   if(i !=3 and i!=5 and i !=0 and i!=1):
 #       m = np.ma.masked_array(m, np.logical_not(mask))
 #       grid_mask = m.mask[grid_pix]
 #       grid_map = np.ma.MaskedArray(m[grid_pix], grid_mask)

    if(i ==4 ):
        m = np.ma.masked_array(m, np.logical_not(map2))
        grid_mask = m.mask[grid_pix]
        grid_map = np.ma.MaskedArray(m[grid_pix], grid_mask)

    if(i ==5 ):
        m = np.ma.masked_array(m, np.logical_not(map4))
        grid_mask = m.mask[grid_pix]
        grid_map = np.ma.MaskedArray(m[grid_pix], grid_mask)


    if(i == 2 or i == 3): vmin,vmax=[-2.,2.]
    if(i == 1): vmin,vmax=[-80,80]
    if(i == 0): vmin,vmax=[-2000,2000]
    if(i == 4 or i == 5): vmin,vmax=[0,600]

    ax = fig.add_subplot(nrow,ncol,i+1, projection='mollweide')
    if(i <=3): 
        image = plt.pcolormesh(longitude[::-1], latitude, m[grid_pix], vmin=vmin,vmax=vmax, rasterized=True, cmap=cmap)
    else:
         image = plt.pcolormesh(longitude[::-1], latitude, m[grid_pix], vmin=vmin,vmax=vmax, rasterized=True, cmap=cmap2)


    ax.set_longitude_grid(60)
    ax.xaxis.set_major_formatter(ThetaFormatterShiftPi(60))
    if width < 22:
        ax.set_latitude_grid(30)
        ax.set_longitude_grid_ends(90)

    #ax.xaxis.set_ticklabels([])
    #ax.yaxis.set_ticklabels([])
    #ax.xaxis.set_ticks([])
    #ax.yaxis.set_ticks([])
    plt.title(title[i], fontsize=12)


    plt.grid()

    if(i <= 3):
        cb = fig.colorbar(image, orientation='horizontal',shrink=.4,
                          pad=0.05, ticks=[vmin, vmax])
        if(i <= 1): cb.ax.xaxis.set_label_text(r'$\mu$K$_{\mathrm {CMB}}$')
        if(i <= 1): cb.ax.xaxis.labelpad = -8
        if(i > 1 and i <=3): cb.ax.xaxis.set_label_text(r'K$_{\mathrm {CMB}}$ rad$^{-2}$')
        if(i > 1 and i <=3): cb.ax.xaxis.labelpad = -8
        # workaround for issue with viewers, see colorbar docstring
        cb.solids.set_edgecolor("face")

    print i
    print vmin,vmax



plt.subplots_adjust(left=0.05, right=0.95, bottom=0.0, top=0.98, wspace=0.10, hspace=0.08)
plt.savefig('pdfs/real_data_6plots.pdf', bbox_inches='tight', pad_inches=cm2inch(0.1))




