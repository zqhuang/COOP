#!/usr/bin/env python
import re
import os
import sys
########by Zhiqi Huang (zqhuang@cita.utoronto.ca) ########
######## This python script hacks cosmomc (any late version >= 2014 April), adding RSD module
#############################################################################


######## If you want to undo the modifications to cosmomc, run the bash script "restore.sh".
#######################################################

coop_propose_updae = 1200
propose_pattern = r'^\s*MPI\_Max\_R\_ProposeUpdate\s*=.*$'
propose_new_pattern = r'^\s*MPI\_Max\_R\_ProposeUpdateNew\s*=.*$'
str_propose = r'MPI_Max_R_ProposeUpdate = '+ str(coop_propose_updae) + r' \nMPI_Max_R_ProposeUpdateNew = '+ str(coop_propose_updae + 200)

def backup_file(fname):
    if(not os.path.isfile(fname + '__.bak')):
        os.system('cp ' + fname + ' ' + fname+'__.bak')


def overwrite_file(fname_from, fname_to):
    if(os.path.isfile(fname_to + '__.bak')):
        os.system('rm -f ' + fname_to)
    else:
        os.system('mv ' + fname_to + ' ' + fname_to + '__.bak')
    os.system('cp ' + fname_from + ' ' + fname_to)

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
    if(sys.argv[1] == "restore"):
        replace_all("source/Makefile", [r'bao\_RSD\.o\s+'], [r'bao.o '])
        os.system('rm -f source/bao_RSD.f90')
        os.system('rm -f test_rsd.ini')
        os.system('rm -f batch2/BAO_RSD.ini')
        sys.exit()

if( not os.path.isfile("source/bao_RSD.f90")):
    print "adding RSD likelihood code"
    if(os.path.isfile("bao_RSD.f90")):
        os.system('cp bao_RSD.f90 source/')
    else:
        print "Cannot find bao_RSD.f90 in the current path"
        print "Aborted."
        sys.exit()
else:
    print "RSD likelihood code already exists"

replace_all("source/Makefile", [r'bao\.o\s+'], [r'bao_RSD.o '])

batch_dir = search_value("test.ini", r'^DEFAULT\((\w+)\/[\w_]*common[\w_]*\.ini\)\s*$')

if batch_dir == '':
    print "Cannot find batch directory."
    print "Aborted"
    sys.exit()


if(os.path.isfile("BAO_RSD.ini")):
    os.system("cp BAO_RSD.ini "+ batch_dir + "/")
else:
    print "Cannot find BAO_RSD.ini"
    print "Aborted"
    sys.exit()


copy_replace_first("test.ini", "test_rsd.ini", [r'^(\#?DEFAULT\(.*BAO.*\))\s*$'], [r'DEFAULT(' + batch_dir + '/BAO_RSD.ini)'])

if(os.path.isfile("BAO_all.ini")):
    os.system("cp BAO_all.ini " + batch_dir + "/")
    copy_replace_first("test.ini", "test_baoall.ini", [r'^(\#?DEFAULT\(.*BAO.*\))\s*$'], [r'DEFAULT(' + batch_dir + '/BAO_all.ini)'])





