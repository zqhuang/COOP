#!/usr/bin/env python
import re
import os
import sys
########by Zhiqi Huang (zqhuang@cita.utoronto.ca) ########
######## This python script hacks cosmomc (any late version > 2013 Oct), adding WL module
################################################################################
def file_is_matched(pattern, fname):
    fp = open(fname, 'r')
    file_content = fp.read()
    fp.close()
    if re.match(r'\^.+\$', pattern) :
        flag = re.I + re.M
    else:
        flag = re.I 
    if re.search(pattern, file_content, flags = flag):
        return True
    else:
        return False


def backup_file(fname):
    if(not os.path.isfile(fname+'__.bak')):
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
    backup_file(fname)
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



def last_line(fname, patterns):
    fp = open(fname, 'r')
    pat = r''
    for line in fp:
        for p in patterns:
            m = re.search(p, line)
            if m:
                pat = p
    return pat

def first_line(fname, patterns):
    fp = open(fname, 'r')
    for line in fp:
        for p in patterns:
            m = re.search(p, line)
            if m:
                return p

def restore(path):
    os.system(r'for i in `ls ' + path + r'/*__.bak`; do cp ${i} ${i/__.bak/}; done')

def abort_quit():
    restore("camb")
    restore("source")
    sys.exit()

#######################################################
if(len(sys.argv)>1):
    if(sys.argv[1] == "restore"  or sys.argv[1] == "remove" or sys.argv[1] == "undo" ):
        replace_all("source/Makefile", [r'\$\(OUTPUT\_DIR\)\/HST\.o\s+\$\(OUTPUT\_DIR\)\/wl\.o'], [r'$(OUTPUT_DIR)/HST.o '])
        replace_all("source/Calculator_CAMB.f90", [r'\ \!\!Zhiqi\ added\[\[[^\]\[]*\!\!\]\]'], [r''])
        replace_all("source/Calculator_Cosmology.f90", [r'\ \!\!Zhiqi\ added\[\[[^\]\[]*\!\!\]\]'], [r''])
        replace_all("source/DataLikelihoods.f90", [r'\ \!\!Zhiqi\ added\[\[[^\]\[]*\!\!\]\]'], [r''])
        replace_first("camb/modules.f90", [r'(integer\s*\,\s*parameter\s*\:\:\s*max\_transfer\_redshifts\s*\=\s*\d+) \+ 100'], [r'\1'])

        os.system('rm -f source/wl.f90')
        os.system('rm -f test_wl.ini')
        os.system('rm -f batch2/WL.ini')
        sys.exit()
    else:
        print sys.argv[1]
        print "Unknown option"
        sys.exit()
else:
    if(os.path.isfile(r"source/wl.f90") and file_is_matched(r'wl\.o', 'source/Makefile')):
        print "wl module is already included."
        print "if you want to remove it, run python link2wl.py restore"
        sys.exit()
            
#first restore to original version
print "adding weak lensing likelihood code"
if(os.path.isfile("wl.f90")):
    if( not os.path.isfile("source/wl.f90") ) :
        os.system("cp wl.f90 source/")
else:
    if( not os.path.isfile("source/wl.f90") ) :
        print "Cannot find wl.f90 in the current path"
        print "Aborted."
        sys.exit()

replace_all("source/Makefile", [r'(\$\(OUTPUT_DIR\)\/HST\.o\s+)'], [r'\1$(OUTPUT_DIR)/wl.o '])


batch_dir = search_value("test.ini", r'^DEFAULT\((\w+)\/[\w_]*common[\w_]*\.ini\)\s*$')

if batch_dir == '':
    print "Cannot find batch directory."
    print "Aborted"
    sys.exit()


if(os.path.isfile("WL.ini")):
    os.system("cp WL.ini "+ batch_dir + "/")
else:
    print "Cannot find WL.ini"
    print "Aborted"
    sys.exit()

copy_replace_first("test.ini", "test_wl.ini", [r'^(\#?DEFAULT\(.*\))\s*$'], [r'\1\nDEFAULT(' + batch_dir + '/WL.ini)'])

replace_first("source/Calculator_CAMB.f90", [r'^(\s*procedure\s*\:\:\s*AngularDiameterDistance\s*\=\>.*)$', r'^(\s*end\s+function\s+CAMBCalc\_AngularDiameterDistance\s*(\!.*)?)$'],[r'\1 !!Zhiqi added[[\n procedure::ComovingRadialDistance => CAMBCalc_ComovingRadialDistance  !!]]', r'\1 !!Zhiqi added[[\n\n   real(mcp) function CAMBCalc_ComovingRadialDistance(this, z) \n     use CAMB, only: ComovingRadialDistance\n     class(CAMB_Calculator) :: this \n     real(mcp),intent(IN):: z \n      CAMBCalc_ComovingRadialDistance = ComovingRadialDistance(z) \n   end function CAMBCalc_ComovingRadialDistance\n!!]]'])


replace_first("source/Calculator_Cosmology.f90", [r'^(\s*procedure\s*\:\:\s*AngularDiameterDistance\s*(\=\>.*)?)$', r'^(\s*end\s+function\s+AngularDiameterDistance\s*(\!.*)?)$'],[r'\1 !!Zhiqi added[[\n procedure::ComovingRadialDistance  !!]]', r'\1 !!Zhiqi added[[\n\n  real(mcp) function ComovingRadialDistance(this, z) \n    class(TCosmologyCalculator) :: this \n     real(mcp),intent(IN):: z \n     call this%ErrorNotImplemented("ComovingRadialDistance")\n        ComovingRadialDistance = 0 \n   end function ComovingRadialDistance\n!!]]'])

replace_first("source/DataLikelihoods.f90", [r'^(\s*subroutine\s+SetDataLikelihoods\s*\(\s*Ini\s*\)\s*(\!.*)?)$', r'^(\s*call\s+BAOLikelihood\_Add\s*\(\s*DataLikelihoods\s*\,\s*Ini\s*\)\s*(!.*)?)$', r'^(\s*CosmoSettings\%use\_LSS\s*\=.*)$'], [r'\1 !!Zhiqi added[[\n   use wl !!]]', r'\1 !!Zhiqi added[[\n   Call WLLikelihood_Add(DataLikelihoods, Ini) !!]]', r'\1 !!Zhiqi added[[\n   CosmoSettings%use_LSS = CosmoSettings%use_LSS .or. use_wl_lss !!]]'])


replace_first("camb/modules.f90", [r'^(\s*integer\s*\,\s*parameter\s*\:\:\s*max\_transfer\_redshifts\s*\=\s*\d+)\s*(\!.*)?$'], [r'\1 + 100'])
