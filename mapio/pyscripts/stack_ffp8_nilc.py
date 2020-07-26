import sys, os, string, math, re

outdir = "tmpmaps/"
spots_dir = "spots/"
stack_dir = "stacked/"
check_files = True

prefix = "ffp8_nilc"
imap_in = "ffp8/" + prefix + "_int_00000_010a_1024.fits"
imask = "planck14/dx11_v2_common_int_mask_010a_1024.fits"
polmap_in = "ffp8/" + prefix + "_pol_00000_hp_20_40_010a_1024.fits"
polmask = "planck14/dx11_v2_common_pol_mask_010a_1024.fits"

fwhm_in = 10
threshold = 0
fwhm_out = 15
unit = "K"  

execfile("stack_common.py")






