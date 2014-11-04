#!/usr/bin/env python
#### Vary data sets and generate cosmomc ini files from a base ini file
#### by  Zhiqi Huang (zqhuang@cita.utoronto.ca)
import re
import os
import sys
import glob
import string


css = [ 'smica', 'commander', 'sevem', 'nilc']
cases = ['case1', 'case2', 'case3', 'case4', 'case4']

for cs in css:
    for case in cases:
        datafile = "planck14/dx11_v2_" + cs + "_pol_" + case + "_cmb_010a_1024.fits"
        ffp8_cmb = "ffp8/dx11_v2_" + cs + "_pol_" + case + "_cmb_mc_00000_010a_1024.fits"
        ffp8_noise = "ffp8/dx11_v2_" + cs + "_pol_" + case + "_noise_mc_00000_010a_1024.fits"
        ffp8file = "ffp8/ffp8_"+cs+"_pol_"+case+"_010a_1024.fits"
        if (os.path.isfile(datafile) and os.path.isfile(ffp8_cmb) and os.path.isfile(ffp8_noise)):
            if(not os.path.isfile(ffp8file)):
                os.system("./MSMAP "+ ffp8_cmb + " " + ffp8_noise + " ADD "+ffp8file)
            os.system("./GetCl " + datafile)
            os.system("./GetCl " + ffp8_cmb)
            os.system("./GetCl " + ffp8_noise)
            os.system("./GetCl " + ffp8file)
            
