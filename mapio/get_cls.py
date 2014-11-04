#!/usr/bin/env python
#### Vary data sets and generate cosmomc ini files from a base ini file
#### by  Zhiqi Huang (zqhuang@cita.utoronto.ca)
import re
import os
import sys
import glob
import string


css = [ 'smica'] #, 'commander', 'sevem', 'nilc']
cases = ['case1'] #, 'case2', 'case3', 'case4', 'case4', 'case5']
mapids = ['00000', '00001', '00002', '00003', '00004',  '00005', '00006', '00007',   '00008',  '00009']

for cs in css:
    for case in cases:
        for mapid in mapids:
            datafile = "planck14/dx11_v2_" + cs + "_pol_" + case + "_cmb_010a_1024.fits"
            ffp8_cmb = "ffp8/dx11_v2_" + cs + "_pol_" + case + "_cmb_mc_" + mapid + "_010a_1024.fits"
            ffp8_noise = "ffp8/dx11_v2_" + cs + "_pol_" + case + "_noise_mc_" + mapid + "_010a_1024.fits"
            ffp8file = "ffp8/ffp8_"+cs+"_pol_"+ case + "_" + mapid + "_010a_1024.fits"
            fig = "clfigs/" + cs + "_" + case + "_mc" + mapid + "_power_comparison.txt"
            if (os.path.isfile(datafile) and os.path.isfile(ffp8_cmb) and os.path.isfile(ffp8_noise)):
                if(not os.path.isfile(ffp8file)):
                    os.system("./MSMAP "+ ffp8_cmb + " " + ffp8_noise + " ADD "+ffp8file)
                    os.system("./GetCl " + datafile)
                    os.system("./GetCl " + ffp8_cmb)
                    os.system("./GetCl " + ffp8_noise)
                    os.system("./GetCl " + ffp8file)
                os.system("./Test " + ffp8_cmb.replace(".fits", "_cls.txt") + " " + ffp8_noise.replace(".fits", "_cls.txt")  + " " + ffp8file.replace(".fits", "_cls.txt")  + " " + datafile.replace(".fits", "_cls.txt")  + " " + fig)
                       # os.system("../utils/fasy.sh " + fig)
            
