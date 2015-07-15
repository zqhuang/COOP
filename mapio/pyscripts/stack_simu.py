import sys, os, string, math, re

outdir = "tmpmaps/"
spots_dir = "spots/"
stack_dir = "stacked/"
check_files = True

prefix = "simu"
imap_in = "simu/simu_int_010a_n1024.fits"
imask = "planck14/dx11_v2_common_int_mask_010a_1024.fits"
polmap_in = "simu/simu_pol_hp_20_40_010a_n1024.fits"
polmask = "planck14/dx11_v2_common_pol_mask_010a_1024.fits"
unit = "K"   

fwhm_in = 10
threshold = 0
fwhm_out = 15

execfile("stack_common.py")




