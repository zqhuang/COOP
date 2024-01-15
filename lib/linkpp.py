#!/usr/bin/env python
import re
import os
import sys
import string
coop_path = r"/mnt/raid-cita/zqhuang/raid-zqhuang/zqhuang/scilibs/COOP"

#######################################################
coop_include = r"-I" + coop_path + r"/include"
coop_link = r"-L" + coop_path + r"/lib" + " -lcoop"
coop_include_append = r"\1 " + coop_include
coop_link_append = r"\1 " + coop_link
coop_propose_updae = 1200
propose_pattern = r'^\s*MPI\_Max\_R\_ProposeUpdate\s*=.*$'
propose_new_pattern = r'^\s*MPI\_Max\_R\_ProposeUpdateNew\s*=.*$'
covmat_pattern = r'^\s*propose\_matrix\s*\=.*$'
covmat_repl = r'propose_matrix = '
str_propose = r'MPI_Max_R_ProposeUpdate = ' + str(coop_propose_updae) + r' \nMPI_Max_R_ProposeUpdateNew = ' + str(coop_propose_updae + 200)

inirootname = r'myinis'


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

def copy_replace_first_append(fname1,  fname2, patterns, repls, app):
    print "copying " + fname1 + " to " + fname2
    fp = open(fname1, 'r')
    file_content = fp.read()
    fp.close()
    for i in range(len(patterns)):
        if(re.match(r'\^.+\$', patterns[i])):
            file_content = re.sub(patterns[i], repls[i], file_content, count = 1,  flags = re.M + re.I)
        else:
            file_content = re.sub(patterns[i], repls[i], file_content, count = 1, flags = re.I)
    if(app != ''):
        file_content = file_content + "\n" + app
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

def subroutine_pattern(func, var):
    return r'subroutine\s+' + func + r'\s*\(\s*' + re.sub(r'\,', r'\s*\,\s*', var) + r'\s*\)\s*(\!.*)?\n(.*\n)+\s*end\s+subroutine\s+' + func


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

def paramnames_list(fname):
    params = []
    fp = open(fname, 'r')
    for line in fp:
        words = string.split(line)
        params.append(words[0])
    fp.close()
    return params

def which_line(params, pattern):
    for i in range(len(params)):
        if(params[i] == pattern):
            return i+1 
    return 0

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

default_params = r'params_CMB.paramnames'
if(not os.path.isfile(default_params)):
    if(os.path.isfile(r'paramnames/params_CMB.paramnames')):
        default_params = r'paramnames/params_CMB.paramnames'
    else:
        print "cannot find the default paramnames file: params_CMB.paramnames"
        sys.exit()
        
param_list = paramnames_list(default_params)

index_logA = which_line(param_list, r'logA')
if index_logA == 0:
    print "Cannot find logA parameter"
    sys.exit()

index_ns = which_line(param_list, r'ns')
if index_ns == 0:
    print "Cannot find ns parameter"
    sys.exit()

index_r = which_line(param_list, r'r')
if index_r == 0:
    print "Cannot find r parameter"
    sys.exit()
    
if numhard + 1 != index_logA:
    print "logA is not the first fast parameter, check the cosmomc version"
    sys.exit()

index_H0 = which_line(param_list, r'H0*')
if index_H0 == 0:
    print "cannot determine the index of H0*"
    sys.exit()


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

if batch_dir == '' or common_file == '':
    print "Cannot find batch directory."
    sys.exit()


print "This version of cosmomc contains:"
print str(index_H0-1) + "theory parameters"
print str(numhard) + " slow parameters"
print str(index_H0-index_logA) + " fast parameters"
print "*****************************************"
print "Modifying files:"

#replace_all("camb/Makefile", [r"^\s*(FFLAGS\s*=.*)$", r"^\s*(F90CRLINK\s*=.*)$"], [ coop_include_append, coop_link_append ])
replace_all("camb/Makefile_main", [r'CAMBOBJ\s*\=\s*\$\(OUTPUT\_DIR\)', r'CAMBSO\s*\=\s*\$\(DLL\_DIR\)'], [r'CAMBOBJ =  $(OUTPUT_DIR)/cambpp.o $(OUTPUT_DIR)',r'CAMBSO =  $(DLL_DIR)/cambpp.o $(DLL_DIR)' ])

#replace_all("source/Makefile", [r"^\s*(FFLAGS\s*=.*)$", r"^\s*(F90CRLINK\s*=.*)$", r"^\s*(INCLUDE\s*=.*)$", r"^\s*(LINKFLAGS\s*=.*)$" ], [ coop_include_append, coop_link_append, coop_include_append, coop_link_append ])

replace_first("camb/power_tilt.f90", [line_pattern(r'module initialpower'), function_pattern('ScalarPower', r'k,ix'), function_pattern('TensorPower', r'k,ix')], [r'module InitialPower\n use camb_mypp', r'function scalarPower(k, ix)\n real(dl) scalarpower, k\n integer ix\n scalarpower = mypp_primordial_ps(k)\n end function scalarpower', r'function tensorPower(k, ix)\n real(dl) tensorPower, k\n integer ix \n tensorPower = mypp_primordial_pt(k) \n end function tensorPower'])

replace_first("source/CosmologyTypes.f90", [line_pattern(r'module cosmologytypes'), line_pattern(r'integer,parameter::max_inipower_params=\d+')], [r'module cosmologyTypes\n use camb_mypp \n', r'integer, parameter:: max_inipower_params = 30'])

replace_first("source/settings.f90", [line_pattern(r'module settings'),  line_pattern(r'integer,parameter::max_theory_params=\d+')], [r'module settings \n use camb_mypp \n', r'integer, parameter:: max_theory_params = 55 \n character(LEN=1024)::cosmomc_paramnames = "'+default_params+'" '])


list_param_pattern = [line_pattern(r'module cosmologyparameterizations'), \
                      r'H0\_min\s*=\s*\d+', \
                      r'H0\_max\s*=\s*\d+', \
                      r"call\s+ini\%read\(\s*\'H0\_min\'\s*\,\s*this\%H0\_min\s*\)", \
                      r"call\s+ini\%read\(\s*\'H0\_max\'\s*\,\s*this\%H0\_max\s*\)", \
                      r'call\s+this\%SetTheoryParameterNumbers\s*\(\s*\d+\,\s*last_power_index\s*\)', \
                      r'^\s*call\s+setfast\s*\(\s*Params\s*\,\s*CMB\s*\)\s*(\!.*)?$', \
                      r"(\'|\")" + re.escape(default_params) + r"(\'|\")", \
                      r'^\s*(real.*\:\:\s*use_min_zre\s*=.*)$', \
                      r'^\s*if\s*\(\s*CMB\%zre\s*\<\s*this\%Use\_min\_zre\)\s*return\s*(\!.*)?$' ]



list_param_replace = [r'module CosmologyParameterizations\n use camb_mypp', \
                      r'H0_min = 55.', \
                      r'H0_max = 85.', \
                      r'this%H0_min = 55.', \
                      r'this%H0_max = 85. \n call Ini%read("use_max_zre", this%use_max_zre)', \
                      r'call this%SetTheoryParameterNumbers('+str(numhard)+', last_power_index + mypp_nknots)', \
                      r'call setfast(params, CMB)\n call mypp_setup_pp(As=exp(params('+str(index_logA)+'))*1.d-10, ns=params('+str(index_ns)+'), r = params('+str(index_r)+'), nknots = mypp_nknots, dlnps = params('+str(numhard)+' + last_power_index+1:'+str(numhard)+' +last_power_index+mypp_nknots))', \
                      r"adjustl(trim(cosmomc_paramnames))", \
                      r'\1\n real(mcp)::use_max_zre = 20.', \
                      r' if (CMB%zre < this%Use_min_zre .or. CMB%zre > this%use_max_zre) return \n ']

replace_all("source/CosmologyParameterizations.f90", list_param_pattern, list_param_replace )

replace_first("source/driver.F90", [line_pattern(r'program cosmomc'), line_pattern(r'call ini%open(inputfile)')], [r'program cosmomc \n \n  use camb_mypp', r'call Ini%Open(InputFile) \n !!start hacking \n mypp_nknots = Ini%Read_Int("mypp_nknots", 0) \n if(mypp_nknots == 0) stop "cannot find mypp_nknots in the ini file." \n cosmomc_paramnames = Ini%Read_String("paramnames", .false.) \n if(trim(cosmomc_paramnames).eq."") stop "cannot find the paramnames entry in the ini file" \n !!end hacking \n '])

print "************************************"
print "Generating params files"


if(not os.path.exists('paramnames')):
    os.system('mkdir paramnames')
if(iniroot != ""):    
    if(not os.path.exists(inirootname)):    
        os.system('mkdir '+ inirootname )

######################   

 
listr = ['100', '050', '010', '001']

for i in range(5, 16, 2):
    ppstr = ''
    for j in range(i):
        ppstr += (r'pp'+ str(j+1) + r'    p_{' + str(j+1) + r'}\n')
    ppstr += r'H0*        H_0'
    copy_replace_first(default_params, r'paramnames/params_scanp' + str(i)+ r'.paramnames',  [ r'^H0\*\s+.+$' ], [ ppstr ])
    copy_replace_first_append(baseini, iniroot + r'dpp'+str(i)+'.ini', [covmat_pattern, common_pattern, r'^file_root\s*=.+$' , r'^action\s*=.+$', r'^compute\_tensors\s*=.+$', r'^(cmb\_dataset\[BICEP2\].*)$', propose_pattern], [covmat_repl, r'DEFAULT(' + batch_dir + '/common_pp.ini) \nmypp_nknots = ' + str(i) + r'\nparamnames = paramnames/params_scanp' + str(i) + r'.paramnames', r'file_root = dpp'+str(i), r'action = 0', r'compute_tensors = T', r'\#\1',str_propose],  r'use_min_zre = 6.' + "\n" + r'use_max_zre = 15.' + "\n\n")
    for fidr in listr:
        copy_replace_first_append(baseini, iniroot + r'dpp'+str(i)+'_fixrp' + fidr + r'.ini', [covmat_pattern, common_pattern, r'^file_root\s*=.+$' , r'^action\s*=.+$', r'^compute\_tensors\s*=.+$', r'^(cmb\_dataset\[BICEP2\].*)$', propose_pattern], [covmat_repl, r'DEFAULT(' + batch_dir + '/common_pp.ini) \npp_model = 1 \nmypp_nknots = ' + str(i) + r'\nparamnames = paramnames/params_scanp' + str(i) + r'.paramnames', r'file_root = dpp'+str(i)+r'_fixrp' + fidr , r'action = 0', r'compute_tensors = T', r'\#\1',str_propose], r'param[r] = 0.' + fidr  + "\n" + r'use_min_zre = 6.' + "\n" + r'use_max_zre = 15.' + "\n\n")

copy_replace_first(common_file, batch_dir + r'/common_pp.ini', [r'^INCLUDE\(params\_CMB\_defaults\.ini\)\s*$'], [r'INCLUDE(params_CMB_pp.ini)'] )
ppstr = r'param[ns] =  0.967  \n'
for i in range(1, 16):
    ppstr += r'param[pp'+ str(i) + r'] = 0. -1. 1. 0.03 0.03 \n'
copy_replace_first_append(batch_dir + r'/params_CMB_defaults.ini', batch_dir + r'/params_CMB_pp.ini', [r'^param\[ns\]\s*=.+$', '^param\[r\]\s*=.+$', '^compute_tensors\s*=.+$'], [ ppstr, r'param[r] = 0.05 0. 0.2 0.03 0.03', r'compute_tensors = T' ], r'use_min_zre = 6.' + "\n" )

