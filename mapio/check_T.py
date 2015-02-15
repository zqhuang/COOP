import time
import os
import sys
cclist = ['commander', 'nilc', 'sevem', 'smica']
nulist = ['0', '1']
hclist = ['hot0', 'cold0', 'hot1', 'cold1']
nslist = ['N', 'S']
readonly = 'F'
nmaps = 0
print "=============== F =============="
for nu in nulist:
    print '------------------------'
    print "nu = " + nu
    for cc in cclist:
        print cc
        for hc in hclist:        
            os.system(r'./SST '+cc+' 1024 T ' + ' ' + nu + ' ' + str(nmaps) + ' self '+ hc + ' T F T ' + readonly + ' > scripts/' + cc + 'T' + nu + 'F.log')
    
sys.exit()
print "=========== NS ==========="

for nu in nulist:
    print '------------------------'    
    print "nu = " + nu
    for cc in cclist:
        print cc
        for ns in nslist:
            os.system(r'./SST ' + cc + ' 1024 T ' + ' ' + nu + ' ' + str(nmaps) + ' self hot0 T  ' + ns + ' T ' + readonly )
            os.system(r'./SST ' + cc + ' 1024 T ' + ' ' + nu + ' ' + str(nmaps) + ' self cold0 T ' + ns + ' T ' + readonly)                      
                    
