#!/usr/bin/env python
#### Vary data sets and generate cosmomc ini files from a base ini file
#### by  Zhiqi Huang (zqhuang@cita.utoronto.ca)
import re
import os
import sys
import glob
import string

inifile = sys.argv[1]
current_path = os.getcwd()
patterns = [r'.*\/', r'scan\_', r'fixrp(\d\d)\d', r'qcdm\_1param', r'qcdm\_3param', 'lowTEB', 'plikTTTEEE', 'plikTT', 'BAO_JLA_HSTlow', 'lens', 'liteTTTEEE', 'liteTT', r'\.ini']
repls = ['', '', r'r\1', 'w1p', 'w3p', 'P', 'E', 'T', 'pr', 'l', 'lE', 'lT', '']
jobname = inifile
for i in range(len(patterns)):
    jobname = re.sub(patterns[i], repls[i], jobname)
    
fp = open(r'scripts/' + jobname + r'.jb', 'w')
fp.write(r'#!/bin/csh -f' + "\n" + r'#PBS -N '+jobname + "\n" + r'#PBS -l nodes=8:ppn=8' + "\n" + r'#PBS -q workq' + "\n" + r'#PBS -l walltime=48:00:00' + "\n" + r'##PBS -r n' + "\n" + r'cd ' + current_path + "\n" + 'mpirun -pernode ./cosmomc ' + inifile + ' > ./scripts/'+jobname+r'.log' + "\n")

fp.close()

os.chdir('scripts')
os.system('qsub ' + jobname + r'.jb')
os.chdir(current_path)


    


