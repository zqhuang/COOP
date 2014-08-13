#!/usr/bin/env python
import re
import os
import sys
########by Zhiqi Huang (zqhuang@cita.utoronto.ca) ########
######## This python script hacks cosmomc (any late version > 2013 Oct), linking it to CosmoRec and adding two parameters A2s1s and tcmb 
#############################################################################

########Before running the script, you need to copy the CosmoRec.v2.0b folder and the patch folder into cosmomc/ directory. Compile CosmoRec first. 
########The script tries to "understand" the cosmomc version you are using. Even if it is slightly modified by you, the script would probably still work. It hacks the source code in both camb/ and source/ directories, then generates two ini files: a2s1s.ini and tcmb.ini. The two ini files are for 6 param + A2s1s and 6 param + T_cmb, respectively. 
####### Re-make cosmomc and submit the ini files. You are done!
###############################################################################


######## If you want to undo the modifications to cosmomc, run the bash script "restore.sh".

patch_path = r'cambcr_patch.CR_v2.0.CAMB_Apr_2014'
#cosmorec_path = r'../CosmoRec.v1.5b'
cosmorec_path = r'CosmoRec.v2.0.1b'
#######################################################

coop_propose_updae = 1200
propose_pattern = r'^\s*MPI\_Max\_R\_ProposeUpdate\s*=.*$'
propose_new_pattern = r'^\s*MPI\_Max\_R\_ProposeUpdateNew\s*=.*$'
str_propose = r'MPI_Max_R_ProposeUpdate = '+ str(coop_propose_updae) + r' \nMPI_Max_R_ProposeUpdateNew = '+ str(coop_propose_updae + 200)

def backup_file(fname):
    if(not os.path.isfile(fname + '__.bak')):
        os.system('cp ' + fname + ' ' + fname+'__.bak')

def replace_first(fname, patterns, repls):
    print "modifying " + fname 
    fp = open(fname, 'r')
    file_content = fp.read()
    fp.close()
    for i in range(len(patterns)):
        if(re.match(r'\^.+\$', patterns[i])):
            file_content = re.sub(patterns[i], repls[i], file_content, count = 1,  flags = re.M + re.I)
        else:
            file_content = re.sub(patterns[i], repls[i], file_content, count = 1, flags = re.I)
    backup_file(fname)
    fp = open(fname, 'w')
    fp.write(file_content)
    fp.close()

def replace_all(fname, patterns, repls):
    print "modifying " + fname 
    fp = open(fname, 'r')
    file_content = fp.read()
    fp.close()
    for i in range(len(patterns)):
        if(re.match(r'\^.+\$', patterns[i])):
            file_content = re.sub(patterns[i], repls[i], file_content, flags = re.M + re.I)
        else:
            file_content = re.sub(patterns[i], repls[i], file_content, flags = re.I)
    backup_file(fname)
    fp = open(fname, 'w')
    fp.write(file_content)
    fp.close()


def copy_replace_first(fname1,  fname2, patterns, repls):
    print "copying " + fname1 + " to " + fname2
    fp = open(fname1, 'r')
    file_content = fp.read()
    fp.close()
    for i in range(len(patterns)):
        if(re.match(r'\^.+\$', patterns[i])):
            file_content = re.sub(patterns[i], repls[i], file_content, count = 1,  flags = re.M + re.I)
        else:
            file_content = re.sub(patterns[i], repls[i], file_content, count = 1, flags = re.I)
    fp = open(fname2, 'w')
    fp.write(file_content)
    fp.close()

def copy_replace_all(fname1, fname2, patterns, repls):
    print "copying " + fname1 + " to " + fname2
    fp = open(fname1, 'r')
    file_content = fp.read()
    fp.close()
    for i in range(len(patterns)):
        if(re.match(r'\^.+\$', patterns[i])):
            file_content = re.sub(patterns[i], repls[i], file_content, flags = re.M + re.I)
        else:
            file_content = re.sub(patterns[i], repls[i], file_content, flags = re.I)
    fp = open(fname2, 'w')
    fp.write(file_content)
    fp.close()


def function_pattern(func, var):
    return r'function\s+' + func + r'\s*\(\s*' + re.sub(r'\,', r'\s*\,\s*', var) + r'\s*\)\s*(\!.*)?\n(.*\n)+\s*end\s+function\s+' + func

def line_pattern(line):
    line = re.sub(r' ', r'\s*', line)
    line = re.sub(r'(\.|\_|\%)', r'\\\1', line)
    line = re.sub(r'\:\:', r'\s*\:\:\s*', line)
    line = re.sub(r'(\=|\,|\(|\))', r'\s*\\\1\s*', line)
    return r'^\s*' + line + '\s*(!.*)?$'

def search_value(fname, pattern):
    fp = open(fname, 'r')
    file_content = fp.read()
    fp.close()
    m = re.search(pattern, file_content, flags = re.M + re.I)
    if m:
        return m.group(1)
    else:
        return ""

def which_line(fname, pattern):
    fp = open(fname, 'r')
    i = 0
    for line in fp:
        i += 1
        m = re.search(pattern, line, flags = re.I)
        if m:
            fp.close()
            return i 
    fp.close()
    return i



def last_line(fname, patterns):
    fp = open(fname, 'r')
    pat = r''
    for line in fp:
        for p in patterns:
            m = re.search(p, line)
            if m:
                pat = p
    return pat

def first_line(fname, patterns):
    fp = open(fname, 'r')
    for line in fp:
        for p in patterns:
            m = re.search(p, line)
            if m:
                return p


def restore(path):
    os.system(r'for i in `ls ' + path + r'/*__.bak`; do cp ${i} ${i/__.bak/}; done')

def abort_quit():
    restore("camb")
    restore("source")
    sys.exit()

#######################################################
#first restore to original version
print "*****************************************"
print "Restoring to original version:"
restore("camb")
restore("source")
if(len(sys.argv) > 1):
    if(sys.argv[1] == "restore"):
        os.system('rm -f camb/*__.bak')
        os.system('rm -f source/*__.bak')
        sys.exit()
print "****************************************"
print "Analyzing the default setups of cosmomc:"
nstr = search_value("source/CosmologyParameterizations.f90",  r'call\s+this\%SetTheoryParameterNumbers\(\s*(\d+)\s*\,\s*last\_power\_index\)') 
if(nstr == ""):
    print "cannot find the number of hard parameters"
    sys.exit()
numhard = int( nstr )

powerpattern = first_line(r"params_CMB.paramnames", [r'^(logA\s+.*)$', r'^(ns\s+.*)$'])

print "Modifying files:"

replace_all(r"source/Makefile", [r"^\s*(RECOMBINATION\s*\?\s*\=).*$", r"^\s*(COSMOREC_PATH\s*\??\=).*$" ], [r'\1cosmorec', r'\1../'+cosmorec_path])

replace_all(r"camb/Makefile_main", [r"^\s*(RECOMBINATION\s*\?\s*\=).*$", r"^\s*(COSMOREC_PATH\s*\??\=).*$" ], [r'\1cosmorec', r'\1../'+cosmorec_path])

backup_file(r'camb/cosmorec.F90')

os.system(r'cp ' + patch_path + r'/cosmorec.F90 camb/cosmorec.F90')

replace_first("source/settings.f90", [r'^\s*module\s+settings\s*(\!.*)?$', line_pattern(r'integer,parameter::max_theory_params=\d+')], [r'module settings\n use constants, only:COBE_CMBTemp', r'integer, parameter:: max_theory_params = 35 \n character(LEN=256)::cosmomc_paramnames = "params_CMB.paramnames" \n integer::cosmomc_num_hard = ' + str(numhard) + r'\n integer::cosmomc_cosmorec_runmode = 0 '] )


replace_first("camb/constants.f90", [r'\,\s*parameter\s*(\:\:\s*COBE_CMBTemp)'], [r'\1'])


replace_first("source/driver.F90", [ line_pattern(r'call ini%open(inputfile)')], [ r'call Ini%Open(InputFile)\n cosmomc_paramnames = Ini%Read_String("paramnames", .false.) \n cosmomc_num_hard = Ini%Read_Int("num_hard", ' + str(numhard) + r') \n cosmomc_cosmorec_runmode = Ini%Read_Int("cosmorec_runmode", 2)'])

replace_all(r"source/Calculator_CAMB.f90", [r'\#ifdef\s+COSMOREC[^\#]*\#else', r'^.*parameter\s*\:\:\s*cons\s*=\s*\(COBE_CMBTemp.*$', r'^(\s*lens\_recon\_scale\s*=\s*CMB\%InitPower.*)$'], [r'#ifdef COSMOREC\n P%Recomb%fdm = CMB%fdm*1.d-23 \n P%Recomb%A2s1s = CMB%A2s1s \n P%Tcmb = COBE_CMBTemp \n P%Recomb%runmode = cosmomc_cosmorec_runmode \n#else', r'real(dl) cons',r'\1 \ncons = (COBE_CMBTemp*1e6)**2'] )

replace_all(r"source/CosmologyTypes.f90", [r'^(\s*Type\s*\,\s*extends.*\:\:\s*CMBParams)\s*(\!.*)?$'], [r'\1\n       real(mcp) A2s1s'])

replace_all(r"source/CosmologyParameterizations.f90",  [r'(call\s+this\%SetTheoryParameterNumbers\(\s*\d+\s*\,\s*last\_power\_index\))', r'(\"|\')params\_CMB\.paramnames(\"|\')', r'^\s*CMB%fdm\s*=\s*Params\((\d+)\)\s*(\!.*)?$'], [r'call this%SetTheoryParameterNumbers(cosmomc_num_hard, last_power_index)', r'trim(cosmomc_paramnames)', r'CMB%fdm = Params(\1) \n CMB%A2s1s = Params(' + str(numhard+1) + r') \n COBE_CMBTemp = Params(' + str(numhard + 2) + ')'])


batch_dir = search_value("test.ini", r'^DEFAULT\((\w+)\/[\w_]*common[\w_]*\.ini\)\s*$')

common_file = search_value("test.ini", r'^DEFAULT\((\w+\/[\w_]*common[\w\_]*\.ini)\)\s*')

common_pattern = r'^(DEFAULT\(\w+\/[\w_]*common[\w\_]*\.ini\))\s*$'


copy_replace_all(r'params_CMB.paramnames', r'params_cosmorec.paramnames', [ powerpattern ], [r'A2s1s        A_{2s\\rightarrow 1s}   #CosmoRec A2s1s parameter \ntcmb        T_{\\rm CMB}   #CosmoRec T_CMB parameter \n\1'] )

copy_replace_first("test.ini", 'a2s1s.ini', [common_pattern, r'^file_root\s*=.+$', r'^action\s*=.+$', propose_pattern], [r'DEFAULT(' + batch_dir + r'/common_a2s1s.ini) \nparamnames = params_cosmorec.paramnames \nnum_hard = ' + str(numhard+2) + r'\ncosmorec_runmode = 0 ', r'file_root = a2s1s', r'action = 0', str_propose] )

copy_replace_all(common_file, batch_dir + r'/common_a2s1s.ini', [r'params\_CMB\_defaults\.ini'],  [r'params_a2s1s.ini'])

copy_replace_all(batch_dir + r'/params_CMB_defaults.ini', batch_dir + r'/params_a2s1s.ini', [r'^param\[fdm\]\s*=.*$'], [r'param[fdm] = 0 0 0 0 0 \nparam[A2s1s] = 8.224 5. 11. 0.8 0.8 \nparam[tcmb] = 2.7255 2.7255 2.7255 0. 0. ' ] )

copy_replace_first("test.ini", 'tcmb.ini', [common_pattern, r'^file_root\s*=.+$', r'^action\s*=.+$', propose_pattern], [r'DEFAULT(' + batch_dir + r'/common_tcmb.ini) \nparamnames = params_cosmorec.paramnames \nnum_hard = '+str(numhard+2) + r'\ncosmorec_runmode = 0 ', r'file_root = tcmb', r'action = 0', str_propose] )

copy_replace_all(common_file, batch_dir + r'/common_tcmb.ini', [r'params\_CMB\_defaults\.ini'],  [r'params_tcmb.ini'])

copy_replace_all(batch_dir + r'/params_CMB_defaults.ini', batch_dir + r'/params_tcmb.ini', [r'^param\[fdm\]\s*=.*$'], [r'param[fdm] = 0 0 0 0 0 \nparam[A2s1s] = 0 0 0 0 0 \nparam[tcmb] = 2.7255 2. 3.5 0.005 0.005' ] )






