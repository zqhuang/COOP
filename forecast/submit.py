#!/usr/bin/env python
### submit cosmomc jobs
#### by  Zhiqi Huang (zqhuang@cita.utoronto.ca)
import re
import os
import sys
import glob
import string

def search_value(fname, pattern):
    fp = open(fname, 'r')
    file_content = fp.read()
    fp.close()
    m = re.search(pattern, file_content, flags = re.M + re.I)
    if m:
        return m.group(1)
    else:
        return ""

if(len(sys.argv) < 2):
    print "Syntax: "
    print "python submit.py IniFile NumChains"
    print "For example:"
    print "python submit.py myinis/qcdm.ini 8"    
    sys.exit()

inifile = sys.argv[1]
if( not os.path.isfile(inifile)):
    print inifile + " does not exist"
    sys.exit()

num_chains = 8    
if(len(sys.argv) >= 3):
    num_chains = int(sys.argv[2])
    if num_chains < 1 or num_chains > 16:
        print "Chain number must be between 1 and 16"
        sys.exit()
    
chainname = search_value(inifile, r'chain\_name\s*\=\s*chains\/([^\s\/\\]+)')
paramnames = search_value(inifile, r'paramnames\s*\=\s*paramnames\/([^\s\/\\]+)')    

if(paramnames ==""):
    print "cannot find paramnames in the ini file " + inifile
    sys.exit()

if(chainname ==""):
    print "cannot find chain_name in the ini file " + inifile
    sys.exit()
    
print "submitting " + inifile
print "chain name: " + chainname
print "paramnames file: " + paramnames

current_path = os.getcwd()

jobname = chainname
    
fp = open(r'scripts/' + jobname + r'.jb', 'w')
fp.write(r'#!/bin/csh -f' + "\n" + r'#PBS -N '+jobname + "\n" + r'#PBS -l nodes=' + str(num_chains/2) + r':ppn=8' + "\n" + r'#PBS -q workq' + "\n" + r'#PBS -l walltime=48:00:00' + "\n" + r'##PBS -r n' + "\n" + r'cd ' + current_path + "\n" + 'mpirun -np ' + str(num_chains) + ' --map-by node:PE=4 ./MCMC ' + inifile + ' > ./scripts/'+jobname+r'.log' + "\n")

fp.close()

os.chdir('scripts')
os.system('qsub ' + jobname + r'.jb')
os.chdir(current_path)


    


