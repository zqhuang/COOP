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
if len(sys.argv) < 3 :
    print "This script scan all the files with the a given extention in a given folder, and determin if the file content matches a given regular expression pattern."
    print "python search.py pattern  file_extention path_to_the_folder groupid"
    print "group_id is only required if you have brackets in the pattern"
    print "for example:"
    print "python search.py '\w+[0-9]{7}' myfile.txt  (search the pattern in myfile.txt in current directory) "
    print "python search.py '\w+[0-9]{7}' myfile.txt ./  (the same as previous one)"  
    print "python search.py 'code(\w+[0-9]{3})' f90 cosmomc/camb 1 (search the pattern in *f90 in cosmomc/camb, print the first group in the bracket"


    sys.exit()
if(force_pattern == ''):
    pattern = sys.argv[1]
else:
    pattern = force_pattern

postfix = sys.argv[2]

if len(sys.argv) < 4 :
    spath = r"./"
else:
    spath = sys.argv[3] + r'/'


if len(sys.argv) < 5 :
    indwant = 0
else:
    indwant = max(0, int(sys.argv[4]))

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
                try:
                    assert isinstance(x, basestring)
                    print x
                except:
                    if( len(x)> indwant):
                        print x[indwant]
                    else:
                        print x[len(x) - 1]

print "*****************************************************************"
print 'match ' + str(nummatch) + ' files in ' + str(nummatch + notmatch) + ' files'
