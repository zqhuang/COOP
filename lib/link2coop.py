#!/usr/bin/env python
import re
import os
import sys
coop_path = r"../../COOP"

#######################################################
coop_include = r"-I" + coop_path + r"/include"
coop_link = r"-L" + coop_path + r"/lib" + " -lcoop"
coop_include_append = r"\1 " + coop_include
coop_link_append = r"\1 " + coop_link

def backup_file(fname):
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

                

#######################################################
#first restore to original version
print "*****************************************"
print "Restoring to original version:"
os.system("./restore.sh")  
print "****************************************"
print "Analyzing the default setups of cosmomc:"
nstr = search_value("source/CosmologyParameterizations.f90",  r'call\s+this\%SetTheoryParameterNumbers\(\s*(\d+)\s*\,\s*last\_power\_index\)') 
if(nstr == ""):
    print "cannot find the number of hard parameters"
    sys.exit()
numhard = int( nstr )

index_w = which_line("params_CMB.paramnames", r'^w\s+.+')
if index_w == 0:
    print "cannot determine the index of dark energy parameter"
    sys.exit()

index_nnu = which_line("params_CMB.paramnames", r'^nnu\s+.+')
if index_nnu == 0:
    print "cannot determine the index of neutrino numbers"
    sys.exit()

index_H0 = which_line("params_CMB.paramnames", r'^H0\*\s+.+')
if index_H0 == 0:
    print "cannot determine the index of H0*"
    sys.exit()


index_logA = which_line("params_CMB.paramnames", r'^logA\s+.+')
if numhard + 1 != index_logA:
    print "Cannot determine the total number of hard parameters"
    sys.exit()


print "This version of cosmomc contains:"
print str(index_H0-1) + "theory parameters"
print str(numhard) + " slow parameters"
print str(index_H0-index_logA) + " fast parameters"
print str(index_nnu-index_w) + " dark energy parameters"
print "*****************************************"
print "Modifying files:"

replace_all("camb/Makefile", [r"^\s*(FFLAGS\s*=.*)$", r"^\s*(F90CRLINK\s*=.*)$"], [ coop_include_append, coop_link_append ]) 

replace_all("source/Makefile", [r"^\s*(FFLAGS\s*=.*)$", r"^\s*(F90CRLINK\s*=.*)$", r"^\s*(INCLUDE\s*=.*)$", r"^\s*(LINKFLAGS\s*=.*)$" ], [ coop_include_append, coop_link_append, coop_include_append, coop_link_append ])


replace_first("camb/equations_ppf.f90", [line_pattern(r"module lambdageneral"), r'is_cosmological_constant\s*=.+', function_pattern('w_de','a'), function_pattern('grho_de', 'a') ], [r'#include "constants.h"\nmodule LambdaGeneral \nuse coop_wrapper', r'is_cosmological_constant = .false.', r'function w_de(a)\n real(dl) w_de, a\n w_de = coop_global_de%wofa(COOP_REAL_OF(a)) \n end function w_de', 'function grho_de(a)\n real(dl) grho_de, a\n grho_de = grhov*coop_global_de%rhoa4_ratio(COOP_REAL_OF(a)) \n end function grho_de'] )

replace_first("camb/power_tilt.f90", [line_pattern(r'module initialpower'), function_pattern('ScalarPower', r'k,ix'), function_pattern('TensorPower', r'k,ix')], [r'module InitialPower\n use coop_wrapper', r'function scalarPower(k, ix)\n real(dl) scalarpower, k\n integer ix\n scalarpower = coop_primordial_ps(k)\n end function scalarpower', r'function tensorPower(k, ix)\n real(dl) tensorPower, k\n integer ix \n tensorPower = coop_primordial_pt(k) \n end function tensorPower'])

replace_first("source/CosmologyTypes.f90", [line_pattern(r'integer,parameter::max_inipower_params=\d+')], [r'integer, parameter:: max_inipower_params = 30'])

replace_first("source/settings.f90", [line_pattern(r'integer,parameter::max_theory_params=\d+')], [r'integer, parameter:: max_theory_params = 50 \n character(LEN=256)::cosmomc_paramnames = "params_CMB.paramnames" '])

replace_all("source/CosmologyParameterizations.f90", [line_pattern(r'module cosmologyparameterizations'), r'H0\_min\s*=\s*\d+', r'H0\_max\s*=\s*\d+', r"call\s+ini\%read\(\s*\'H0\_min\'\s*\,\s*this\%H0\_min\s*\)", r"call\s+ini\%read\(\s*\'H0\_max\'\s*\,\s*this\%H0\_max\s*\)", r'call\s+this\%SetTheoryParameterNumbers\s*\(\s*\d+\,\s*last_power_index\s*\)', r'end\s+subroutine\s+setforH', r'^\s*CMB\%nnu\s*=\s*params\(\d+\)\s*$', r'^\s*CMB\%YHe\s*=\s*params\(\d+\)\s*$', r'^\s*CMB\%iso\_cdm\_correlated\s*=\s*params\(\d+\)\s*$', r'^\s*CMB\%zre\_delta\s*=\s*params\(\d+\)\s*$', r'^\s*CMB\%Alens\s*=\s*params\(\d+\)\s*$',  r'^\s*CMB\%Alensf\s*=\s*params\(\d+\)\s*$',  r'^\s*CMB\%fdm\s*=\s*params\(\d+\)\s*$', r'^\s*call\s+setfast\s*\(\s*Params\s*\,\s*CMB\s*\)\s*$', r"'params_CMB.paramnames'"], [r'#include "constants.h"\nmodule CosmologyParameterizations\n use coop_wrapper', r'H0_min = 55.', r'H0_max = 85.', r'this%H0_min = 55.', r'this%H0_max = 85.', r'call this%SetTheoryParameterNumbers(cosmomc_de_index + cosmomc_de_num_params+cosmomc_de2pp_num_params-1, cosmomc_pp_num_params)', r'call coop_setup_cosmology_from_cosmomc(params, H0/100.d0)\nend subroutine setForH', r'CMB%nnu = Params(cosmomc_de_index + cosmomc_de_num_params)', r'CMB%YHe = Params(cosmomc_de_index + cosmomc_de_num_params+1)',  r'CMB%iso_cdm_correlated = Params(cosmomc_de_index + cosmomc_de_num_params+2)', r'CMB%zre_delta = Params(cosmomc_de_index + cosmomc_de_num_params+3)', r'CMB%Alens = Params(cosmomc_de_index + cosmomc_de_num_params+4)', r'CMB%Alensf = Params(cosmomc_de_index + cosmomc_de_num_params+5)', r'CMB%fdm = Params(cosmomc_de_index + cosmomc_de_num_params+6)', r"call setfast(params, CMB)\n call coop_setup_cosmology_from_cosmomc(params)\n call coop_setup_pp()", r"adjustl(trim(cosmomc_paramnames))" ])

replace_first("source/driver.F90", [line_pattern(r'program cosmomc'), line_pattern(r'call ini%open(inputfile)')], [r'program cosmomc \n  use coop_wrapper', r'call Ini%Open(InputFile)\n cosmomc_de_model = Ini%Read_Int("de_model", cosmomc_de_model) \n cosmomc_de_num_params = Ini%Read_Int("de_num_params", cosmomc_de_num_params)\n cosmomc_pp_model = Ini%Read_Int("pp_model", cosmomc_pp_model) \n cosmomc_pp_num_params = Ini%Read_Int("pp_num_params", cosmomc_pp_num_params) \n cosmomc_paramnames = Ini%Read_String("paramnames", .false.) \n if(trim(cosmomc_paramnames).eq."") stop "cannot find the paramnames entry in the ini file" \n cosmomc_de_index = ' + str(index_w) + r'\n cosmomc_de2pp_num_params = ' + str(numhard - 1) + r' - cosmomc_de_index \n coop_pp_cosmomc_num = last_power_index'])

print "************************************"
print "Generating params files"

batch_dir = search_value("test.ini", r'^DEFAULT\((\w+)\/[\w_]*common[\w_]*\.ini\)\s*$')
common_file = search_value("test.ini", r'^DEFAULT\((\w+\/[\w_]*common[\w\_]*\.ini)\)\s*')
common_pattern = r'^(DEFAULT\(\w+\/[\w_]*common[\w\_]*\.ini\))\s*$'

if batch_dir == '' or common_file == '':
    print "Cannot find batch directory."
    sys.exit()

os.system('mkdir paramnames')
os.system('cp params_CMB.paramnames paramnames/')
copy_replace_first("params_CMB.paramnames", "paramnames/params_qcdm.paramnames", [r'^w\s+.+$', r'^wa\s+.+$'], [r"epss          \epsilon_s", r"epsinf           \epsilon_{\infty}\nzetas           \zeta_s" ] )


##lcdm
copy_replace_first("test.ini", 'lcdm.ini', [common_pattern, r'^file_root\s*=.+$', r'^action\s*=.+$'], [r'\1 \nde_model = 2 \nde_num_params='+str(index_nnu-index_w)+r'\npp_model = 0 \npp_num_params = ' + str(index_H0 - index_logA) +r'\nparamnames = paramnames/params_CMB.paramnames', r'file_root = lcdm', r'action = 0'] )


### w0
copy_replace_first("test.ini", 'w0.ini', [common_pattern,  r'^file_root\s*=.+$', r'^action\s*=.+$'], [r'DEFAULT(' + batch_dir + r'/common_w0.ini) \nde_model = 2 \nde_num_params='+str(index_nnu-index_w)+r'\npp_model = 0 \npp_num_params = ' + str(index_H0 - index_logA) +r'\nparamnames = paramnames/params_CMB.paramnames', r'file_root = w0', r'action = 0'] )
copy_replace_first(common_file, batch_dir + r'/common_w0.ini', [r'^INCLUDE\(params_CMB_defaults\.ini\)\s*$'], [r'INCLUDE(params_CMB_w0.ini)'] )
copy_replace_first(batch_dir + r'/params_CMB_defaults.ini', batch_dir + r'/params_CMB_w0.ini', [r'^param\[w\]\s*=.+$', r'^param\[wa\]\s*=.+$'], [ r'param[w] = -1 -3 1 0.05 0.05', r'param[wa] = 0 0 0 0 0' ] ) 

copy_replace_first("test.ini", 'w0wa.ini', [common_pattern, r'^file_root\s*=.+$', r'^action\s*=.+$'], [r'DEFAULT(' + batch_dir + r'/common_w0wa.ini) \nde_model = 2 \nde_num_params='+str(index_nnu-index_w)+r'\npp_model = 0 \npp_num_params = ' + str(index_H0 - index_logA) +r'\nparamnames = paramnames/params_CMB.paramnames', r'file_root = w0wa', r'action = 0'] )
copy_replace_first(common_file, batch_dir + r'/common_w0wa.ini', [r'^INCLUDE\(params_CMB_defaults\.ini\)\s*$'], [r'INCLUDE(params_CMB_w0wa.ini)'] )
copy_replace_first(batch_dir + r'/params_CMB_defaults.ini', batch_dir + r'/params_CMB_w0wa.ini', [r'^param\[w\]\s*=.+$', r'^param\[wa\]\s*=.+$'], [ r'param[w] = -1 -3 1 0.05 0.05', r'param[wa] = 0. -3 3 0.2 0.2' ] ) 


#qcdm 1 parameter: epss
copy_replace_first('test.ini', 'epss.ini', [common_pattern, r'^file_root\s*=.+$', r'^action\s*=.+$'], [r'DEFAULT(' + batch_dir + r'/common_epss.ini) \nde_model = 3\nde_num_params=3\nparamnames = paramnames/params_qcdm.paramnames', r'file_root = epss', r'action = 0'] )
copy_replace_first(common_file, batch_dir + r'/common_epss.ini', [r'^INCLUDE\(params_CMB_defaults\.ini\)\s*$'], [r'INCLUDE(params_CMB_epss.ini)'] )
copy_replace_first(batch_dir + r'/params_CMB_defaults.ini', batch_dir + r'/params_CMB_epss.ini', [r'^param\[w\]\s*=.+$'], [ r'param[w] = -1 -1 -1 0 0 \nparam[epss] = 0 -1.5 1.5 0.1 0.1 \nparam[epsinf] = 0 0 0 0 0  \nparam[zetas] = 0 0 0 0 0' ] ) 


##qcdm
copy_replace_first('test.ini', 'qcdm.ini', [common_pattern, r'^file_root\s*=.+$', r'^action\s*=.+$'], [r'DEFAULT('  + batch_dir + r'/common_qcdm.ini) \nde_model = 3\nde_num_params=3\nparamnames = paramnames/params_qcdm.paramnames', r'file_root = qcdm', r'action = 0'] )
copy_replace_first(common_file, batch_dir + r'/common_qcdm.ini', [r'^INCLUDE\(params_CMB_defaults\.ini\)\s*$'], [r'INCLUDE(params_CMB_qcdm.ini)'] )
copy_replace_first(batch_dir + r'/params_CMB_defaults.ini', batch_dir + r'/params_CMB_qcdm.ini', [r'^param\[w\]\s*=.+$'], [ r'param[w] = -1 -1 -1 0 0 \nparam[epss] = 0 -1.5 1.5 0.1 0.1 \nparam[epsinf] = 0.05 0 1. 0.05 0.05 \nparam[zetas] = 0 -1 1 0.1 0.1' ] ) 
 
ppnum = index_H0 - index_logA
for i in range(7, 16):
    ppstr = ''
    for j in range(i):
        ppstr += (r'pp'+ str(j+1) + r'    p_{' + str(j+1) + r'}\n')
    ppstr += r'H0*        H_0'
    copy_replace_first(r'params_CMB.paramnames', r'paramnames/params_scanp' + str(i)+ r'.paramnames',  [ r'^H0\*\s+.+$' ], [ ppstr ])
    copy_replace_first('test.ini', 'scanp'+str(i)+'.ini', [common_pattern, r'^file_root\s*=.+$' , r'^action\s*=.+$', r'^compute\_tensors\s*=.+$', r'^\#(cmb\_dataset\[BICEP2\].*)$'], [r'DEFAULT(' + batch_dir + '/common_pp.ini) \npp_model = 1 \npp_num_params = ' + str(ppnum+i) + r'\nparamnames = paramnames/params_scanp' + str(i) + r'.paramnames', r'file_root = pp'+str(i) , r'action = 0', r'compute_tensors = T', r'\1'])


copy_replace_first(common_file, batch_dir + r'/common_pp.ini', [r'^INCLUDE\(params_CMB_defaults\.ini\)\s*$'], [r'INCLUDE(params_CMB_pp.ini)'] )
ppstr = r'param[ns] =  0.96 0.96 0.96 0 0 \n'
for i in range(1, 16):
    ppstr += r'param[pp'+ str(i) + r'] = 0. -1. 1. 0.03 0.03 \n'
copy_replace_first(batch_dir + r'/params_CMB_defaults.ini', batch_dir + r'/params_CMB_pp.ini', [r'^param\[ns\]\s*=.+$', '^param\[r\]\s*=.+$', '^compute_tensors\s*=.+$'], [ ppstr, r'param[r] = 0.12 0. 1. 0.05 0.05', r'compute_tensors = T' ] ) 





