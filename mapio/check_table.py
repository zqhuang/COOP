import time
import os
cclist = ['commander', 'nilc', 'sevem', 'smica']
nulist = ['0', '1']
hlist = ['hot0', 'hot1', 'hot2']
clist = ['cold0', 'cold1', 'cold2']
nslist = ['N', 'S']
readonly = 'T'
nmaps = 300
print "=============== F =============="
for nu in nulist:
    print " ------------------------- "    
    print "nu = " + nu
    print " * hot spots * "
    for cc in cclist:
        print cc
        for h in hlist:        
            os.system(r'./SST '+ cc + r' 1024 QU ' + nu + r' ' + str(nmaps) + r' self '+ h + r' T F T ' + readonly)
    print(" * cold spots * ")
    for cc in cclist:
        print cc
        for c in clist:        
            os.system(r'./SST '+ cc + r' 1024 QU ' + nu + r' ' + str(nmaps) + r' self '+ c + r' T F T '+ readonly)


print "=============== NS =============="

for nu in nulist:
    print " ------------------------- "
    print "nu = " + nu
    for cc in cclist:
        print cc
        for ns in nslist:
            print "hemisphere: " + ns
            os.system(r'./SST ' + cc + ' 1024 QU ' + ' ' + nu + ' 0 self hot0 T  ' + ns + ' T '+ readonly)
            os.system(r'./SST ' + cc + ' 1024 QU ' + ' ' + nu + ' 0 self cold0 T ' + ns + ' T '+ readonly)                      
                    
