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


print "Modifying files:"

replace_all(r"source/Makefile", [r"^\s*(RECOMBINATION\s*\?\s*\=).*$", r"^\s*(COSMOREC_PATH\s*\??\=).*$" ], [r'\1cosmorec', r'\1../'+cosmorec_path])

replace_all(r"camb/Makefile_main", [r"^\s*(RECOMBINATION\s*\?\s*\=).*$", r"^\s*(COSMOREC_PATH\s*\??\=).*$" ], [r'\1cosmorec', r'\1../'+cosmorec_path])

backup_file(r'camb/cosmorec.F90')

os.system(r'cp ' + patch_path + r'/cosmorec.F90 camb/cosmorec.F90')

replace_all(r"source/Calculator_CAMB.f90", [r'\#ifdef\s+COSMOREC[^\#]*\#else'], [r'#ifdef COSMOREC\n P%Recomb%fdm = CMB%fdm*1.d-23 \n P%Recomb%A2s1s = CMB%A2s1s \n if(P%Recomb%fdm .gt. 0.) P%Recomb%runmode = 3 \n#else'] )

replace_all(r"source/CosmologyTypes.f90", [r'^(\s*Type\s*\,\s*extends.*\:\:\s*CMBParams)\s*$'], [r'\1\n       real(mcp) A2s1s'])

replace_all(r"source/CosmologyParameterizations.f90",  [r'(call\s+this\%SetTheoryParameterNumbers\(\s*\d+\s*\,\s*last\_power\_index\))', r'params\_CMB\.paramnames', r'^\s*CMB%fdm\s*=\s*Params\((\d+)\)\s*$'], [r'call this%SetTheoryParameterNumbers(' + str(numhard + 1) + r', last_power_index)', r'params_cosmorec.paramnames', r'CMB%fdm = Params(\1) \n CMB%A2s1s = Params(\1 + 1)'])

copy_replace_all(r'params_CMB.paramnames', r'params_cosmorec.paramnames', [r'^(fdm\s+.*)$'], [r'\1 \nA2s1s        A_{2s\\rightarrow 1s}   #CosmoRec A2s1s parameter'] )




batch_dir = search_value("test.ini", r'^DEFAULT\((\w+)\/[\w_]*common[\w_]*\.ini\)\s*$')

common_file = search_value("test.ini", r'^DEFAULT\((\w+\/[\w_]*common[\w\_]*\.ini)\)\s*')

common_pattern = r'^(DEFAULT\(\w+\/[\w_]*common[\w\_]*\.ini\))\s*$'

copy_replace_all(batch_dir + r'/params_CMB_defaults.ini', batch_dir + r'/params_cosmorec_defaults.ini', [r'^param\[fdm\]\s*=.*$'], [r'param[fdm] = 0 0 0 0 0 \nparam[A2s1s] = 8.224 7 10 0.01 0.01' ] )

copy_replace_all(common_file, batch_dir + r'/common_cosmorec.ini', [r'params\_CMB\_defaults\.ini'],  [r'params_cosmorec_defaults.ini'])


copy_replace_first("test.ini", 'a2s1s.ini', [common_pattern, r'^file_root\s*=.+$', r'^action\s*=.+$'], [r'DEFAULT(' + batch_dir + r'/common_cosmorec.ini) \nparamnames = params_cosmorec.paramnames', r'file_root = a2s1s', r'action = 0'] )

