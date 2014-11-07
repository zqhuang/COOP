#!/usr/bin/env python
#### Vary data sets and generate cosmomc ini files from a base ini file
#### by  Zhiqi Huang (zqhuang@cita.utoronto.ca)
import re
import os
import sys
import glob
import string

cs = 'smica'
mapids = ['00000', '00001', '00002', '00003', '00004',  '00005', '00006', '00007',   '00008',  '00009']

for mapid in mapids:
    cmb_imap = "ffp8/dx11_v2_" + cs + "_int_cmb_mc_" + mapid + "_010a_1024.fits"
    cmb_polmap = "ffp8/dx11_v2_" + cs + "_pol_case1_cmb_mc_" + mapid + "_hp_20_40_010a_1024.fits"
    noise_imap = "ffp8/dx11_v2_" + cs + "_int_noise_mc_" + mapid + "_010a_1024.fits"
    noise_polmap = "ffp8/dx11_v2_" + cs + "_pol_case1_noise_mc_" + mapid + "_hp_20_40_010a_1024.fits"
    cmb_output = "ffp8/ffp8_iqu_cmb_" + mapid + "_cls.txt"
    noise_output = "ffp8/ffp8_iqu_noise_" + mapid + "_cls.txt"
    os.system("./GetCl " + cmb_imap + " " + cmb_polmap+ " " +  cmb_output)
    os.system("./GetCl " + noise_imap + " " + noise_polmap+ " " +  noise_output)
    
            
