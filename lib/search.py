#!/usr/bin/env python
import re
import os
import sys
import glob

force_pattern = ''
########by Zhiqi Huang (zqhuang@cita.utoronto.ca) ########
def file_match(pattern, fname):
    fp = open(fname, 'r')
    file_content = fp.read()
    fp.close()
    if re.match(r'\^.+\$', pattern) :
        flag = re.I + re.M
    else:
        flag = re.I 
    return re.findall(pattern, file_content, flags = flag)

        
#######################################################
if len(sys.argv) <= 1 :
    print "This script scan all the files with the a given extention in a given folder, and determin if the file content matches a given regular expression pattern."
    print "python search.py pattern  file_extention path_to_the_folder"
    print "for example:"
    print "python search.py '\w+[0-9]{3}' f90 cosmomc/camb"
    sys.exit()
if(force_pattern == ''):
    pattern = sys.argv[1]
else:
    pattern = force_pattern

if len(sys.argv) <= 2:
    postfix=r'*'
else:
    postfix = sys.argv[2]

if len(sys.argv) <= 3 :
    spath = r"./"
else:
    spath = sys.argv[3] + r'/'

if(os.path.isfile(spath+postfix)):
    prefix=''
else:
    prefix=r'*.'

print 'pattern:'
print pattern
nummatch = 0
notmatch = 0
print "********************* match list *******************************"
for fname in glob.glob(spath + prefix  + postfix):
    if os.path.isfile(fname):
        res =  file_match(pattern, fname)
        if len(res) == 0:
            notmatch += 1
        else:
            print fname 
            nummatch += 1
            for x in res:
                print "\t" + x
if(postfix !=  postfix.upper()):
    for fname in glob.glob(spath + prefix + postfix.upper()):
        if os.path.isfile(fname):
            res =  file_match(pattern, fname)
            if len(res) == 0:
                notmatch += 1
            else:
                print fname 
                nummatch += 1
                for x in res:
                    print "\t" + x

print "*****************************************************************"
print 'match ' + str(nummatch) + ' files in ' + str(nummatch + notmatch) + ' files'
