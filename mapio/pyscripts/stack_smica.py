import sys, os, string, math, re

outdir = "tmpmaps/"
spots_dir = "peaks/"
stack_dir = "stacked/"
check_files = True

prefix = "smica"

fwhm_in = 10
nside = 2048 * 5 / fwhm_in
if(nside < 1000):
    strnside = "0"+str(nside)
else:
    strnside = str(nside)
    

imap_in = "planck14/dx11_v2_"+prefix+"_int_cmb_0" + str(fwhm_in) + "a_" + strnside + ".fits"
imask = "planck14/dx11_v2_common_int_mask_0" + str(fwhm_in) + "a_" + strnside + ".fits"
polmap_in = "planck14/dx11_v2_"+prefix+"_pol_case1_cmb_hp_20_40_0" + str(fwhm_in) + "a_" + strnside + ".fits"
polmask = "planck14/dx11_v2_common_pol_mask_0" + str(fwhm_in) + "a_" + strnside + ".fits"

threshold = 0
fwhm_out = 15
#fwhm_in

execfile("stack_common.py")





