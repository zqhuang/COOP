#!/usr/bin/env python
import re
import os
import sys
import glob

force_pattern = ''
########by Zhiqi Huang (zqhuang@cita.utoronto.ca) ########
#######################################################


if len(sys.argv) < 4 :
    print "This script scan all the files with the a given extention in a given folder, replace all the occrence of the pattern with replace_string."
    print "python replace.py pattern  replace_string file_extention [path_to_the_folder]"
    print "for example:"
    print "python replace.py '\w+[0-9]{7}' 'xxx' myfile.txt  (replace the pattern with xxx in *myfile.txt in current directory) "
    print "python search.py '\w+[0-9]{7}' 'xxx' myfile.txt ./  (the same as previous one)"  
    print "python search.py 'code(\w+[0-9]{3})' 'xxx' f90 cosmomc/camb  (replace the pattern with xxx in *f90 in cosmomc/camb, print the first group in the bracket"
    sys.exit()


        
if(force_pattern == ''):
    pattern = sys.argv[1]
else:
    pattern = force_pattern
    
repl_str = sys.argv[2]
repl_str = repl_str.replace(r'\n', "\n")
postfix = sys.argv[3]

if len(sys.argv) <5 :
    spath = r"./"
else:
    spath = sys.argv[4] + r'/'


if(os.path.isfile(spath+postfix)):
    prefix=''
else:
    prefix=r'*.'

print 'pattern:'
print pattern

def do_repl(mobj):
    return repl_str;

def file_replace(pattern, fname):
    fp = open(fname, 'r')
    file_content = fp.read()
    fp.close()
    if re.match(r'\^.+\$', pattern) :
        flag = re.I + re.M
    else:
        flag = re.I    
    new_content = re.sub(pattern, do_repl, file_content, flags = flag)
    if(new_content != file_content):
        print fname + " is modified"
        os.system("cp " + fname + " " + fname + "__.bak")
        fp = open(fname, 'w')
        fp.write(new_content)
        fp.close()



for fname in glob.glob(spath + prefix  + postfix):
    if os.path.isfile(fname):
        file_replace(pattern, fname)

