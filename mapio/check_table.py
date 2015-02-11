import time
import os
cclist = ['commander', 'nilc', 'sevem', 'smica']
nulist = ['0', '1']
hlist = ['hot0', 'hot1', 'hot2']
clist = ['cold0', 'cold1', 'cold2']
nslist = ['N', 'S']
print "=============== F =============="
for nu in nulist:
    print "nu = " + nu
    print " hot spots"
    for cc in cclist:
        print cc
        for h in hlist:        
            os.system(r'./SST '+cc+' 1024 QU ' + ' ' + nu + ' 0 self '+ h + ' T F T')
    print(" cold spots")
    for cc in cclist:
        print cc
        for c in clist:        
            os.system(r'./SST '+cc+' 1024 QU ' + ' ' + nu + ' 0 self '+ c + ' T F T')


print "=========== NS ==========="

for nu in nulist:
    print "nu = " + nu
    for cc in cclist:
        print cc
        for ns in nslist:
            print "hemisphere: " + ns
            os.system(r'./SST ' + cc + ' 1024 QU ' + ' ' + nu + ' 0 self hot0 T  ' + ns + ' T')
            os.system(r'./SST ' + cc + ' 1024 QU ' + ' ' + nu + ' 0 self cold0 T ' + ns + ' T')                      
                    
