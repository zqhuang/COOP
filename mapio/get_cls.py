#!/usr/bin/env python
#### Vary data sets and generate cosmomc ini files from a base ini file
#### by  Zhiqi Huang (zqhuang@cita.utoronto.ca)
import re
import os
import sys
import glob
import string

postfix = "_hp_20_40_010a_1024.fits"

css = [ 'smica'] #, 'commander', 'sevem', 'nilc']
mapids = ['00000', '00001', '00002', '00003', '00004',  '00005', '00006', '00007',   '00008',  '00009']

for cs in css:
    for mapid in mapids:
        datafile = "planck14/dx11_v2_" + cs + "_pol_case1_cmb" + postfix
        ffp8_cmb = "ffp8/dx11_v2_" + cs + "_pol_case1_cmb_mc_" + mapid + "_010a_1024.fits"
        ffp8_noise = "ffp8/dx11_v2_" + cs + "_pol_case1_noise_mc_" + mapid + "_010a_1024.fits"
        ffp8file = "ffp8/ffp8_pol_" + mapid + "_010a_1024.fits"
        fig = "clfigs/" + cs + "_mc" + mapid + "_power_comparison.txt"
        if (os.path.isfile(datafile) and os.path.isfile(ffp8_cmb) and os.path.isfile(ffp8_noise)):
            if(not os.path.isfile(ffp8file)):
                os.system("./MSMAP "+ ffp8_cmb + " " + ffp8_noise + " ADD "+ffp8file)
                os.system("./GetCl " + datafile)
                os.system("./GetCl " + ffp8_cmb)
                os.system("./GetCl " + ffp8_noise)
                os.system("./GetCl " + ffp8file)
            os.system("./FAB " + ffp8_cmb.replace(".fits", "_cls.txt") + " " + ffp8_noise.replace(".fits", "_cls.txt")  + " " + ffp8file.replace(".fits", "_cls.txt")  + " " + datafile.replace(".fits", "_cls.txt")  + " " + fig)
            # os.system("../utils/fasy.sh " + fig)
            
