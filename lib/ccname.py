#!/usr/bin/env python
import re
import os
import sys
import string

numchains = 8
plist = [r".converge_stat", r".inputparams", r".paramnames", r".ranges", r".likelihoods"]

namefrom = "chains/" + sys.argv[1]
nameto =  "chains/" + sys.argv[2]

def rename_file(ffrom, fto):
    if(not os.path.isfile(ffrom)):
        print ffrom + " is missing"
        return
    if(os.path.isfile(fto)):
        print  fto + " already exists"
        return
    os.system("mv " + ffrom + " " + fto)
    return


for postfix in plist:
    rename_file(namefrom + postfix, nameto + postfix)

for i in range(1, numchains+1):
    rename_file(namefrom + "_" + str(i) + ".txt", nameto + "_" + str(i) + ".txt")
    rename_file(namefrom + "_" + str(i) + ".log", nameto + "_" + str(i) + ".log")
    rename_file(namefrom + "_" + str(i) + ".chk", nameto + "_" + str(i) + ".chk")

