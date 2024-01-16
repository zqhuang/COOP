import sys, os, string, math, re

outdir = "tmpmaps/"
spots_dir = "ffp8_spots/"
stack_dir = "ffp8_stacked/"
check_files = True

prefix = "smica"
mapid  = "00000"

fwhm_in = 10
nside = 2048 * 5 / fwhm_in
if(nside < 1000):
    strnside = "0"+str(nside)
else:
    strnside = str(nside)
    

imap_in = "ffp8/dx11_v2_"+prefix+"_int_cmb_mc_"+mapid+"_0" + str(fwhm_in) + "a_" + strnside + ".fits"
imask = "planck14/dx11_v2_common_int_mask_0" + str(fwhm_in) + "a_" + strnside + ".fits"
polmap_in = "ffp8/dx11_v2_"+prefix+"_pol_case1_cmb_mc_"+mapid+"_hp_20_40_0" + str(fwhm_in) + "a_" + strnside + ".fits"
polmask = "planck14/dx11_v2_common_pol_mask_0" + str(fwhm_in) + "a_" + strnside + ".fits"

threshold = 0
fwhm_out = 15

execfile("stack_common.py")






