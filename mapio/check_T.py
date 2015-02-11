import time
import os
cclist = ['commander', 'nilc', 'sevem', 'smica']
nulist = ['0', '1']
hclist = ['hot0', 'cold0', 'hot1', 'cold1']
nslist = ['N', 'S']
print "=============== F =============="
for nu in nulist:
    print "nu = " + nu
    print " hot spots"
    for cc in cclist:
        print cc
        for hc in hclist:        
            os.system(r'./SST '+cc+' 1024 T ' + ' ' + nu + ' 0 self '+ hc + ' T F T')
    print(" cold spots")
    

print "=========== NS ==========="

for nu in nulist:
    print "nu = " + nu
    for cc in cclist:
        print cc
        for ns in nslist:
            print "hemisphere: " + ns
            os.system(r'./SST ' + cc + ' 1024 T ' + ' ' + nu + ' 0 self hot0 T  ' + ns + ' T')
            os.system(r'./SST ' + cc + ' 1024 T ' + ' ' + nu + ' 0 self cold0 T ' + ns + ' T')                      
                    
