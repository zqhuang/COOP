#!/usr/bin/env python
import re
import os
import sys
patch_path = r'cambcr_patch.CR_v2.0.CAMB_Apr_2014'
#cosmorec_path = r'../CosmoRec.v1.5b'
cosmorec_path = r'CosmoRec.v2.0b'

#######################################################

def backup_file(fname):
    os.system('cp ' + fname + ' ' + fname+'__.bak')

def replace_first(fname, patterns, repls):
    print "modifying " + fname 
    fp = open(fname, 'r')
    file_content = fp.read()
    fp.close()
    for i in range(len(patterns)):
        if(re.match(r'\^.+\$', patterns[i])):
            file_content = re.sub(patterns[i], repls[i], file_content, count = 1,  flags = re.M + re.I)
        else:
            file_content = re.sub(patterns[i], repls[i], file_content, count = 1, flags = re.I)
    backup_file(fname)
    fp = open(fname, 'w')
    fp.write(file_content)
    fp.close()

def replace_all(fname, patterns, repls):
    print "modifying " + fname 
    fp = open(fname, 'r')
    file_content = fp.read()
    fp.close()
    for i in range(len(patterns)):
        if(re.match(r'\^.+\$', patterns[i])):
            file_content = re.sub(patterns[i], repls[i], file_content, flags = re.M + re.I)
        else:
            file_content = re.sub(patterns[i], repls[i], file_content, flags = re.I)
    os.system('cp ' + fname + ' ' + fname+'__.bak')
    fp = open(fname, 'w')
    fp.write(file_content)
    fp.close()


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

def copy_replace_all(fname1, fname2, patterns, repls):
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
    fp.write(file_content)
    fp.close()


def function_pattern(func, var):
    return r'function\s+' + func + r'\s*\(\s*' + re.sub(r'\,', r'\s*\,\s*', var) + r'\s*\)\s*(\!.*)?\n(.*\n)+\s*end\s+function\s+' + func

def line_pattern(line):
    line = re.sub(r' ', r'\s*', line)
    line = re.sub(r'(\.|\_|\%)', r'\\\1', line)
    line = re.sub(r'\:\:', r'\s*\:\:\s*', line)
    line = re.sub(r'(\=|\,|\(|\))', r'\s*\\\1\s*', line)
    return r'^\s*' + line + '\s*(!.*)?$'

def search_value(fname, pattern):
    fp = open(fname, 'r')
    file_content = fp.read()
    fp.close()
    m = re.search(pattern, file_content, flags = re.M + re.I)
    if m:
        return m.group(1)
    else:
        return ""

def which_line(fname, pattern):
    fp = open(fname, 'r')
    i = 0
    for line in fp:
        i += 1
        m = re.search(pattern, line, flags = re.I)
        if m:
            fp.close()
            return i 
    fp.close()
    return i

#######################################################
#first restore to original version
print "*****************************************"
print "Restoring to original version:"
os.system("./restore.sh")  
print "****************************************"
print "Analyzing the default setups of cosmomc:"
nstr = search_value("source/CosmologyParameterizations.f90",  r'call\s+this\%SetTheoryParameterNumbers\(\s*(\d+)\s*\,\s*last\_power\_index\)') 
if(nstr == ""):
    print "cannot find the number of hard parameters"
    sys.exit()
numhard = int( nstr )



index_w = which_line("params_CMB.paramnames", r'^w\s+.+')
if index_w == 0:
    print "cannot determine the index of dark energy parameter"
    sys.exit()

index_nnu = which_line("params_CMB.paramnames", r'^nnu\s+.+')
if index_nnu == 0:
    print "cannot determine the index of neutrino numbers"
    sys.exit()

index_H0 = which_line("params_CMB.paramnames", r'^H0\*\s+.+')
if index_H0 == 0:
    print "cannot determine the index of H0*"
    sys.exit()


index_logA = min(which_line("params_CMB.paramnames", r'^logA\s+.+'), which_line("params_CMB.paramnames", r'^ns\s+.+'))

if numhard + 1 != index_logA:
    print "Cannot determine the total number of hard parameters"
    sys.exit()


print "This version of cosmomc contains:"
print str(index_H0-1) + "theory parameters"
print str(numhard) + " slow parameters"
print str(index_H0-index_logA) + " fast parameters"
print str(index_nnu-index_w) + " dark energy parameters"
print "*****************************************"
print "Modifying files:"

replace_all("source/Makefile", [r"^\s*RECOMBINATION\s*\??\=.*$" ], [r'RECOMBINATION=cosmorec'])

replace_all("camb/Makefile_main", [r"^\s*RECOMBINATION\s*\??\=.*$", r"^\s*COSMOREC_PATH\s*\??\=.*$" ], [r'RECOMBINATION=cosmorec', r'COSMOREC_PATH=../'+cosmorec_path])

backup_file(r'camb/cosmorec.F90')
os.system(r'cp ' + patch_path + r'/cosmorec.F90 camb/cosmorec.F90')





