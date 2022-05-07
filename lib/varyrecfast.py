#!/usr/bin/env python
import re
import os
import sys
########by Zhiqi Huang (zqhuang@cita.utoronto.ca) ########
######## This python script hacks cosmomc (any late version > 2013 Oct), adding two parameters A2s1s and tcmb. A2s1s is used in Recfast, and tcmb is used globally (including C_l normalization).
#############################################################################
do_cl_norm = False

coop_propose_updae = 1200
propose_pattern = r'^\s*MPI\_Max\_R\_ProposeUpdate\s*=.*$'
propose_new_pattern = r'^\s*MPI\_Max\_R\_ProposeUpdateNew\s*=.*$'
str_propose = r'MPI_Max_R_ProposeUpdate = '+ str(coop_propose_updae) + r' \nMPI_Max_R_ProposeUpdateNew = '+ str(coop_propose_updae + 200)

inirootname = "myinis"

if(inirootname != ""):
    iniroot = inirootname + r'/'
else:
    iniroot = ""


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

if(not os.path.exists('paramnames')):
    os.system('mkdir paramnames')
default_params = r'params_CMB.paramnames'
if(not os.path.isfile(default_params)):
    if(os.path.isfile(r'paramnames/params_CMB.paramnames')):
        default_params = r'paramnames/params_CMB.paramnames'
    else:
        print "cannot find the default paramnames file: params_CMB.paramnames"
        sys.exit()

    
if(iniroot != ""):    
    if(not os.path.exists(inirootname)):    
        os.system('mkdir '+ inirootname )


powerpattern = first_line(default_params, [r'^(logA\s+.*)$', r'^(ns\s+.*)$'])

print "Modifying files:"


replace_first("source/settings.f90", [r'^\s*module\s+settings\s*(\!.*)?$', line_pattern(r'integer,parameter::max_theory_params=\d+')], [r'module settings\n use constants, only:COBE_CMBTemp', r'integer, parameter:: max_theory_params = 35 \n character(LEN=256)::cosmomc_paramnames = "'+default_params+r'" \n integer::cosmomc_num_hard = ' + str(numhard+2) ] )


replace_first("camb/constants.f90", [r'\,\s*parameter\s*(\:\:\s*COBE_CMBTemp)'], [r'\1'])


replace_first("source/driver.F90", [ line_pattern(r'call ini%open(inputfile)')], [ r'call Ini%Open(InputFile)\n cosmomc_paramnames = Ini%Read_String("paramnames", .false.) \n cosmomc_num_hard = Ini%Read_Int("num_hard", ' + str(numhard+2) + r') '])

if(do_cl_norm):
    replace_all(r"source/Calculator_CAMB.f90", [r'^(\s*subroutine\s+CAMBCalc_CMBToCAMB\s*\(.+\).*)$', r'^(\s*P\%YHe\s*\=\s*CMB\%YHe.*)$', r'^.*parameter\s*\:\:\s*cons\s*=\s*\(COBE_CMBTemp.*$', r'^(\s*lens\_recon\_scale\s*=\s*CMB\%InitPower.*)$'], [r'\1 \n use recdata, only:Lambda', r'\1 \n Lambda = CMB%A2s1s \n P%Tcmb = COBE_CMBTemp', r'real(dl) cons',r'\1 \ncons = (COBE_CMBTemp*1e6)**2'] )
else:
    replace_all(r"source/Calculator_CAMB.f90", [r'^(\s*subroutine\s+CAMBCalc_CMBToCAMB\s*\(.+\).*)$', r'^(\s*P\%YHe\s*\=\s*CMB\%YHe.*)$', r'^(.*parameter\s*\:\:\s*cons\s*=\s*)\(COBE_CMBTemp.*$'], [r'\1 \n use recdata, only:Lambda', r'\1 \n Lambda = CMB%A2s1s \n P%Tcmb = COBE_CMBTemp', r'\1 (2.72558d6)**2'] )

replace_all(r"source/CosmologyTypes.f90", [r'^(\s*Type\s*\,\s*extends.*\:\:\s*CMBParams)\s*(\!.*)?$'], [r'\1\n       real(mcp) A2s1s'])

replace_all(r"source/CosmologyParameterizations.f90",  [r'(call\s+this\%SetTheoryParameterNumbers\(\s*\d+\s*\,\s*last\_power\_index\))', r'(\"|\')'+re.escape(default_params)+'(\"|\')', r'^\s*CMB%fdm\s*=\s*Params\((\d+)\)\s*(\!.*)?$'], [r'call this%SetTheoryParameterNumbers(cosmomc_num_hard, last_power_index)', r'trim(adjustl(cosmomc_paramnames))', r'CMB%fdm = Params(\1) \n CMB%A2s1s = Params(' + str(numhard+1) + r') \n COBE_CMBTemp = Params(' + str(numhard + 2) + ')'])

baseini = 'my_base.ini'

if(not os.path.isfile(baseini)):
    if(os.path.isfile("test.ini")):
        patterns = [r'^\s*checkpoint\s*\=.*$']
        repls = [ r'checkpoint = T']
        copy_replace_first('test.ini',  baseini, patterns, repls)
    else:
        print "test.ini file does not exist."
        sys.exit()


batch_dir = search_value(baseini, r'^DEFAULT\((\w+)\/[\w_]*common[\w_]*\.ini\)\s*$')

common_file = search_value(baseini, r'^DEFAULT\((\w+\/[\w_]*common[\w\_]*\.ini)\)\s*')

common_pattern = r'^(DEFAULT\(\w+\/[\w_]*common[\w\_]*\.ini\))\s*$'


copy_replace_all(default_params, r'paramnames/params_recfast.paramnames', [ powerpattern ], [r'A2s1s        A_{2s\\rightarrow 1s}   # A2s1s parameter \ntcmb        T_0   # T_CMB parameter \n\1'] )


copy_replace_first(baseini, iniroot + r'recfast_a2s1s.ini', [r'^propose\_matrix\s*\=.*$', common_pattern, r'^file_root\s*=.+$', r'^action\s*=.+$', propose_pattern], [r'propose_matrix = plots/recfast_a2s1s.covmat', r'DEFAULT(' + batch_dir + r'/common_a2s1s.ini) \nparamnames = paramnames/params_recfast.paramnames \nnum_hard = ' + str(numhard+2) + r'\nrecfast_runmode = 0 ', r'file_root = recfast_a2s1s', r'action = 0', str_propose] )


copy_replace_all(common_file, batch_dir + r'/common_a2s1s.ini', [r'params\_CMB\_defaults\.ini'],  [r'params_a2s1s.ini'])

copy_replace_all(batch_dir + r'/params_CMB_defaults.ini', batch_dir + r'/params_a2s1s.ini', [r'^param\[fdm\]\s*=.*$'], [r'param[fdm] = 0  \nparam[A2s1s] = 8.224 4. 11. 0.6 0.6 \nparam[tcmb] = 2.7255 ' ] )

copy_replace_first(baseini, iniroot + 'recfast_tcmb.ini', [r'^propose\_matrix\s*\=.*$', common_pattern, r'^file_root\s*=.+$', r'^action\s*=.+$', propose_pattern], [r'propose_matrix = plots/recfast_tcmb.covmat', r'DEFAULT(' + batch_dir + r'/common_tcmb.ini) \nparamnames = paramnames/params_recfast.paramnames \nnum_hard = '+str(numhard+2) + r'\nrecfast_runmode = 0 ', r'file_root = recfast_tcmb', r'action = 0', str_propose] )

copy_replace_all(common_file, batch_dir + r'/common_tcmb.ini', [r'params\_CMB\_defaults\.ini'],  [r'params_tcmb.ini'])

copy_replace_all(batch_dir + r'/params_CMB_defaults.ini', batch_dir + r'/params_tcmb.ini', [r'^param\[fdm\]\s*=.*$'], [r'param[fdm] = 0  \nparam[A2s1s] = 8.2245809 \nparam[tcmb] = 2.72558 1.5 5. 0.01 0.01' ] )


copy_replace_first(baseini, iniroot + r'recfast_nnu.ini', [r'^propose\_matrix\s*\=.*$', common_pattern, r'^file_root\s*=.+$', r'^action\s*=.+$', propose_pattern], [r'propose_matrix = plots/recfast_nnu.covmat', r'DEFAULT(' + batch_dir + r'/common_nnu.ini) \nparamnames = paramnames/params_recfast.paramnames \nnum_hard = '+str(numhard+2) + r'\nrecfast_runmode = 0 ', r'file_root = recfast_nnu', r'action = 0', str_propose] )

copy_replace_all(common_file, batch_dir + r'/common_nnu.ini', [r'params\_CMB\_defaults\.ini'],  [r'params_nnu.ini'])

copy_replace_all(batch_dir + r'/params_CMB_defaults.ini', batch_dir + r'/params_nnu.ini', [r'^param\[fdm\]\s*=.*$', r'^param\[nnu\]\s*=.*$'], [r'param[fdm] = 0  \nparam[A2s1s] = 8.2245809 \nparam[tcmb] = 2.72558', r'param[nnu] = 3.046 1.1 5.1 0.3 0.3'] )

copy_replace_first(baseini, iniroot+r'recfast_yhe.ini', [r'^propose\_matrix\s*\=.*$', common_pattern, r'^file_root\s*=.+$', r'^action\s*=.+$', propose_pattern], [r'propose_matrix = plots/recfast_yhe.covmat', r'DEFAULT(' + batch_dir + r'/common_yhe.ini) \nparamnames = paramnames/params_recfast.paramnames \nnum_hard = '+str(numhard+2) + r'\nrecfast_runmode = 0 ', r'file_root = recfast_yhe', r'action = 0', str_propose] )

copy_replace_all(common_file, batch_dir + r'/common_yhe.ini', [r'params\_CMB\_defaults\.ini'],  [r'params_yhe.ini'])

copy_replace_all(batch_dir + r'/params_CMB_defaults.ini', batch_dir + r'/params_yhe.ini', [r'^param\[fdm\]\s*=.*$', r'^param\[yhe\]\s*=.*$', r'^bbn_consistency\s*\=.*$'], [r'param[fdm] = 0  \nparam[A2s1s] = 8.2245809 \nparam[tcmb] = 2.72558', r'param[yhe] =  0.25 0.2 0.3 0.01 0.01', r'bbn_consistency = F'] )


copy_replace_all(batch_dir + r'/params_CMB_defaults.ini', batch_dir + r'/params_fdm.ini', [r'^param\[fdm\]\s*=.*$'], [r'param[fdm] = 0.05 0. 1. 0.05 0.05  \nparam[A2s1s] = 8.2245809 \nparam[tcmb] = 2.72558'] )


copy_replace_first(baseini, iniroot + r'recfast_lcdm.ini', [ common_pattern, r'^file_root\s*=.+$', r'^action\s*=.+$', propose_pattern], [ r'DEFAULT(' + batch_dir + r'/common_lcdm.ini) \nparamnames = paramnames/params_recfast.paramnames \nnum_hard = ' + str(numhard+2) + r'\nrecfast_runmode = 0 ', r'file_root = recfast_lcdm', r'action = 0', str_propose] )

copy_replace_all(common_file, batch_dir + r'/common_lcdm.ini', [r'params\_CMB\_defaults\.ini'],  [r'params_lcdm.ini'])

copy_replace_all(batch_dir + r'/params_CMB_defaults.ini', batch_dir + r'/params_lcdm.ini', [r'^param\[fdm\]\s*=.*$'], [r'param[fdm] = 0  \nparam[A2s1s] = 8.224  \nparam[tcmb] = 2.72558 ' ] )


