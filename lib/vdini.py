#!/usr/bin/env python
#### Vary data sets and generate cosmomc ini files from a base ini file
#### by  Zhiqi Huang (zqhuang@cita.utoronto.ca)
import re
import os
import sys
import glob
import string


###  define the short names for data files
datafile = dict()
datafile['lowTEB'] = r'batch2/lowTEB.ini'
datafile['plikTT'] = r'batch2/plik_dx11dr2_HM_v16_TT.ini'
datafile['plikTTTEEE'] = r'batch2/plik_dx11dr2_HM_v16_TTTEEE.ini'
datafile['lens'] = r'batch2/lensing.ini'
datafile['BAO'] = r'batch2/BAO.ini'
datafile['JLA'] = r'batch2/JLA.ini'
datafile['HSTlow'] = r'batch2/HST_GPE70p6.ini'
datafile['HSThigh'] = r'batch2/HST_high.ini'
datafile['WL'] = r'batch2/WL.ini'
datafile['RSD'] =  r'batch2/BAO_RSD.ini'




###  Edit the collection of data sets that you want to run
want = [ ['lowTEB', 'plikTT'], \
         ['lowTEB', 'plikTT', 'BAO', 'JLA', 'HSTlow'], \
         ['lowTEB', 'plikTT', 'RSD'], \
         ['lowTEB', 'plikTT', 'WL'], \
         ['lowTEB', 'plikTT', 'RSD', 'WL'], \
         ['lowTEB', 'plikTT', 'BAO', 'JLA', 'HSThigh'], \
         ['lowTEB', 'plikTT', 'lens'], \
         ['lowTEB', 'plikTT', 'lens', 'BAO', 'JLA', 'HSTlow'], \
         ['lowTEB', 'plikTT', 'lens', 'RSD'], \
         ['lowTEB', 'plikTT', 'lens', 'WL'], \
         ['lowTEB', 'plikTT', 'lens', 'RSD', 'WL'], \
         ['lowTEB', 'plikTT', 'lens', 'BAO', 'JLA', 'HSThigh'], \
         ['lowTEB', 'plikTTTEEE'], \
         ['lowTEB', 'plikTTTEEE', 'BAO', 'JLA', 'HSTlow'] ]

         

def copy_replace_first(fname1,  fname2, patterns, repls):
    print "copying " + fname1 + " to " + fname2
    fp = open(fname1, 'r')
    file_content = fp.read()
    fp.close()
    for i in range(len(patterns)):
        if(re.match(r'\^.+\$', patterns[i])):
            file_content = re.sub(patterns[i], repls[i], file_content, count = 1,  flags = re.M + re.I)
        else:
            file_content = re.sub(patterns[i], repls[i], file_content, count = 1, flags = re.I)
    fp = open(fname2, 'w')
    fp.write(file_content)
    fp.close()

def copy_replace_all(fname1, fname2, patterns, repls, headstr):
    print "copying " + fname1 + " to " + fname2
    fp = open(fname1, 'r')
    file_content = fp.read()
    fp.close()
    for i in range(len(patterns)):
        if(re.match(r'\^.+\$', patterns[i])):
            file_content = re.sub(patterns[i], repls[i], file_content, flags = re.M + re.I)
        else:
            file_content = re.sub(patterns[i], repls[i], file_content, flags = re.I)
    fp = open(fname2, 'w')
    fp.write(headstr + "\n" + file_content)
    fp.close()


def generate_ini(fname, datasets):
    postfix = ''
    inc = ''
    for data in datasets:
        postfix = postfix + "_" + data
        inc = inc + r'DEFAULT('  + datafile[data] + r')' + "\n"
    fout = fname.replace('.ini', postfix + '.ini')
    patterns = [r'^\s*(DEFAULT\([^\(\)]*\))\s*$', \
                r'^\#\#\_OUT(DEFAULT\(batch[^\(\)]*common[^\(\)]*\.ini\))$', \
                r'^\s*file_root\s*\=([\w\d\_\-\.]*)\s*$']
    repls = [r'##_OUT\1', \
             r'\1', \
             r'file_root = \1'+postfix ]
    copy_replace_all(fname, fout, patterns, repls, inc)

fname = sys.argv[1]
if (not os.path.isfile(fname)):
    print fname + " does not exist"
    sys.exit()

    
for datasets in want:
    generate_ini(fname, datasets)


    


