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

default_covmat = r'planck_covmats/base_planck_lowl_lowLike.covmat'  ##this covmat will be used by default

################# No need to change anything below ####################

def copy_replace_all(fname1, fname2, patterns, repls, headstr):
    print "copying " + fname1 + " to " + fname2
    fp = open(fname1, 'r')
    file_content = fp.read()
    fp.close()
    isdist = re.search(r'plot\_data\_dir',file_content, flags=re.I + re.M)
    if(isdist is None):
        pstart = 0
    else:
        pstart = 2
    for i in range(pstart, len(patterns)):
        if(re.match(r'\^.+\$', patterns[i])):
            file_content = re.sub(patterns[i], repls[i], file_content, flags = re.M + re.I)
        else:     
            file_content = re.sub(patterns[i], repls[i], file_content, flags = re.I)
    fp = open(fname2, 'w')
    if(pstart == 0 ):
        fp.write(headstr + "\n" + file_content)
    else:
        fp.write(file_content)
    fp.close()

def search_value(fname, pattern):
    fp = open(fname, 'r')
    file_content = fp.read()
    fp.close()
    m = re.search(pattern, file_content, flags = re.M + re.I)
    if m:
        return m.group(1)
    else:
        return ""


def generate_ini(fname, datasets, action):
    postfix = ''
    inc = ''
    for data in datasets:
        postfix = postfix + "_" + data
        inc = inc + r'DEFAULT('  + datafile[data] + r')' + "\n"
    fout = fname.replace('.ini', postfix + '.ini')
    patterns = [ r'^\s*(DEFAULT\([^\(\)]*\))\s*$', \
                r'^\#\#\_OUT(DEFAULT\(batch[^\(\)]*common[^\(\)]*\.ini\))$', \
                 r'^\s*file\_root\s*\=\s*(\S*)\s*$']
    repls = [r'##_OUT\1', \
             r'\1', \
             r'file_root = \1'+postfix ]
    froot = search_value(fname, r'^\s*file\_root\s*=\s*(\S*)\s*$')
    covf = search_value(fname, r'^\s*propose\_matrix\s*=\s*' + r'((\S*(\/|\\))?' + froot + r')\.covmat\s*$')
    if(covf != ''):
        newcovf = covf + postfix + r'.covmat'
        if(os.path.isfile(newcovf)):
            patterns.append(r'^\s*propose\_matrix\s*=\s*\S*\s*$')
            repls.append(r'propose_matrix = ' + newcovf)
        elif(os.path.isfile(covf)):
            patterns.append(r'^\s*propose\_matrix\s*=\s*\S*\s*$')
            repls.append(r'propose_matrix = ' + covf)
        elif(os.path.isfile(default_covmat)):
            patterns.append(r'^\s*propose\_matrix\s*=\s*\S*\s*$')
            repls.append(r'propose_matrix = ' + default_covmat)
        else:
            patterns.append(r'^\s*propose\_matrix\s*=\s*\S*\s*$')
            repls.append(r'propose_matrix = ')
    copy_replace_all(fname, fout, patterns, repls, inc)
    if(action != ""):
        os.system(action +" " + fout)

if(len(sys.argv) < 2):
    print "syntax:"
    print  "python vdini.py  YourIniFileName"
    sys.exit()

    
fname = sys.argv[1]

if(len(sys.argv)>=3):
    action = sys.argv[2]
else:
    action = ""

if (not os.path.isfile(fname)):
    print fname + " does not exist"
    sys.exit()

    
for datasets in want:
    generate_ini(fname, datasets, action)

    


