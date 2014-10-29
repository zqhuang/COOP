#!/usr/bin/env python
import re
import os
import sys
import string

numchains = 8
plist = [r".converge_stat", r".inputparams", r".paramnames", r".ranges", r".likelihoods"]

namefrom = "chains/" + sys.argv[1]
nameto =  "chains/" + sys.argv[2]


if os.path.isfile(nameto + ".inputparams"):
    print "the target name already exists"
    sys.exit()

for postfix in plist:
    if(os.path.isfile( namefrom + postfix)):
        os.system("mv " + namefrom + postfix + " " + nameto + postfix)

for i in range(1, numchains+1):
    os.system("mv " + namefrom + "_" + str(i) + ".txt " + nameto + "_" + str(i) + ".txt")
    os.system("mv " + namefrom + "_" + str(i) + ".log " + nameto + "_" + str(i) + ".log")
    if(os.path.isfile(namefrom + "_" + str(i) + ".chk ")):
        os.system("mv " + namefrom + "_" + str(i) + ".chk " + nameto + "_" + str(i) + ".chk")
