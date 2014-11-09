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
covmat_repl = r'propose_matrix = plots/'
str_propose = r'MPI_Max_R_ProposeUpdate = '+ str(coop_propose_updae) + r' \nMPI_Max_R_ProposeUpdateNew = '+ str(coop_propose_updae + 200)

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

if(os.path.isfile('params_CMB.paramnames')):
    param_list = paramnames_list("params_CMB.paramnames")
else:
    print "cannot find the default paramnames file: params_CMB.paramnames"
    sys.exit()

index_w = which_line(param_list, r'w')
if index_w == 0:
    print "cannot determine the index of dark energy parameter"
    sys.exit()

index_wa = which_line(param_list, r'wa')
if index_wa == 0:
    num_w_params = 1
else:
    num_w_params = 2

index_post_w = index_w + num_w_params

index_logA = which_line(param_list, r'logA')
if index_logA == 0:
    print "Cannot find logA parameter"
    sys.exit()

if numhard + 1 != index_logA:
    print "logA is not the first fast parameter"
    print "This is not compatible with the default settings in COOP"
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
print str(num_w_params) + " dark energy parameters"
print "*****************************************"
print "Modifying files:"

replace_all("camb/Makefile", [r"^\s*(FFLAGS\s*=.*)$", r"^\s*(F90CRLINK\s*=.*)$"], [ coop_include_append, coop_link_append ]) 

replace_all("source/Makefile", [r"^\s*(FFLAGS\s*=.*)$", r"^\s*(F90CRLINK\s*=.*)$", r"^\s*(INCLUDE\s*=.*)$", r"^\s*(LINKFLAGS\s*=.*)$" ], [ coop_include_append, coop_link_append, coop_include_append, coop_link_append ])

replace_all("camb/modules.f90", [r'^(\s*subroutine\s+inithermo\s*\(\s*taumin\s*\,\s*taumax\s*\)\s*(\!.*)?)$', r'call\s+Recombination\_Init\s*\(\s*CP\%Recomb\s*\,\s*CP\%omegac\s*\,\s*CP\%omegab'], [r'\1\n use coop_wrapper', 'call Recombination_Init(CP%Recomb, CP%Omegac*coop_global_cdm%density_ratio(1.d0/1090.d0)/1090.d0**3, CP%Omegab*coop_global_cdm%density_ratio(1.d0/1090.d0)/1090.d0**3'])

replace_all("camb/equations_ppf.f90", \
            [ line_pattern(r"module lambdageneral"),  \
              r'is_cosmological_constant\s*=.+', function_pattern('w_de','a'), \
              function_pattern('grho_de', 'a'), line_pattern(r"module GaugeInterface"), \
              r'\=\s*grhoc\s*\/\s*a\s*\n', \
              r'\=\s*grhob\s*\/\s*a\s*\n' , \
              r'\(\s*grhoc\s*\+\s*grhob\s*\)\s*\*\s*a\s*\+', \
              r'\(\s*grhoc\s*\+\s*grhob\s*\)\s*\/\s*a\s*\n', \
              r'\(\s*clxc\s*\*\s*grhoc\s*\+\s*clxb\s*\*\s*grhob\s*\)\s*\/\s*a\s*\n', \
              r'\(\s*grhob\s*\+\s*grhoc\)\/sqrt', \
              r'^\s*call\s+Nu_rho\s*\(\s*a\s*\*\s*nu\_masses\s*\(\s*nu\_i\s*\)\s*\,\s*rhonu\)\s*(\!.*)?$' , \
              r'^\s*rhonudot\s*\=\s*Nu_drho\s*\(\s*a\s*\*\s*nu\_masses\s*\(\s*nu\_i\s*\)\s*\,\s*adotoa\s*\,\s*rhonu\)\s*(\!.*)?$', \
              r'^\s*call\s+Nu\_background\s*\(\s*a\s*\*\s*nu\_masses\s*\(\s*nu\_i\s*\)\s*\,\s*rhonu\s*\,\s*pnu\s*\)\s*(\!.*)?$'], \
            [ r'#include "constants.h"\nmodule LambdaGeneral \nuse coop_wrapper', \
              r'is_cosmological_constant = .false.', \
              r'function w_de(a)\n real(dl) w_de, a\n w_de = coop_global_de%wofa(COOP_REAL_OF(a)) \n end function w_de', \
              r'function grho_de(a)\n real(dl) grho_de, a\n grho_de = grhov*coop_global_de%rhoa4_ratio(COOP_REAL_OF(a)) \n end function grho_de', \
              r'module GaugeInterface \n use coop_wrapper', \
              r'   = grhoc * coop_global_cdm%density_ratio(a)*a**2 \n ', \
              r'   = grhob * coop_global_baryon%density_ratio(a)*a**2 \n ', \
              r'(grhoc*coop_global_cdm%rhoa4_ratio(a)+grhob*coop_global_baryon%rhoa4_ratio(a)) + ', \
              r'(grhoc*coop_global_cdm%density_ratio(a)+grhob*coop_global_baryon%density_ratio(a))*a**2 \n', \
              r'(clxc*grhoc*coop_global_cdm%density_ratio(a) + clxb*grhob*coop_global_baryon%density_ratio(a))*a**2 \n', \
              r'(grhob + grhoc)*(coop_global_cdm%density_ratio(coop_min_scale_factor)*coop_min_scale_factor**3)/sqrt', \
              r'rhonu = coop_global_massive_neutrinos%Omega/coop_global_massive_neutrinos%Omega_massless * coop_global_massive_neutrinos%rhoa4_ratio(a)', \
              r'rhonudot = -3.d0*coop_global_massive_neutrinos%wp1effofa(a) * rhonu * adotoa', \
              r'    rhonu = coop_global_massive_neutrinos%Omega/coop_global_massive_neutrinos%Omega_massless * coop_global_massive_neutrinos%rhoa4_ratio(a) \n    pnu = rhonu*coop_global_massive_neutrinos%wofa(a)'] )

replace_first("camb/power_tilt.f90", [line_pattern(r'module initialpower'), function_pattern('ScalarPower', r'k,ix'), function_pattern('TensorPower', r'k,ix')], [r'module InitialPower\n use coop_wrapper', r'function scalarPower(k, ix)\n real(dl) scalarpower, k\n integer ix\n scalarpower = coop_primordial_ps(k)\n end function scalarpower', r'function tensorPower(k, ix)\n real(dl) tensorPower, k\n integer ix \n tensorPower = coop_primordial_pt(k) \n end function tensorPower'])

replace_first("source/CosmologyTypes.f90", [line_pattern(r'integer,parameter::max_inipower_params=\d+')], [r'integer, parameter:: max_inipower_params = 30'])

replace_first("source/settings.f90", [line_pattern(r'integer,parameter::max_theory_params=\d+')], [r'integer, parameter:: max_theory_params = 50 \n character(LEN=256)::cosmomc_paramnames = "params_CMB.paramnames" '])

replace_all("camb/halofit_ppf.f90",  \
            [ r'^(\s*subroutine\s+NonLinear_GetNonLinRatios\s*\(\s*CAMB\_Pk\)\s*(\!.*)?)$',  \
              r'^\s*om\_m\s*\=\s*omega\_m\s*\(.*$', r'^\s*om\_v\s*\=\s*omega\_v\s*\(.*$', \
              r'\(\s*1\.?\s*\+\s*w\_lam\s*\)' ],  \
            [ r'\1\n use coop_wrapper', \
              r'  om_m = (coop_global_cdm%Omega * coop_global_cdm%rhoa4_ratio(a) + coop_global_baryon%Omega * coop_global_baryon%rhoa4_ratio(a))/coop_global_cosmology%H2a4(a)',  \
              r' om_v = coop_global_de%Omega * coop_global_de%rhoa4_ratio(a) / coop_global_cosmology%H2a4(a) ', \
              r' coop_global_de%wp1 ' ])

list_param_pattern = [line_pattern(r'module cosmologyparameterizations'), \
                      r'H0\_min\s*=\s*\d+', \
                      r'H0\_max\s*=\s*\d+', \
                      r"call\s+ini\%read\(\s*\'H0\_min\'\s*\,\s*this\%H0\_min\s*\)", \
                      r"call\s+ini\%read\(\s*\'H0\_max\'\s*\,\s*this\%H0\_max\s*\)", \
                      r'call\s+this\%SetTheoryParameterNumbers\s*\(\s*\d+\,\s*last_power_index\s*\)', \
                      r'end\s+subroutine\s+setforH',  \
                      r'^\s*call\s+setfast\s*\(\s*Params\s*\,\s*CMB\s*\)\s*(\!.*)?$', \
                      r"(\'|\")params_CMB.paramnames(\'|\")"]



list_param_replace = [r'#include "constants.h"\nmodule CosmologyParameterizations\n use coop_wrapper', \
                      r'H0_min = 55.', \
                      r'H0_max = 85.', \
                      r'this%H0_min = 55.', \
                      r'this%H0_max = 85.', \
                      r'call this%SetTheoryParameterNumbers(cosmomc_de_index + cosmomc_de_num_params+cosmomc_de2pp_num_params-1, cosmomc_pp_num_params)', \
                      r'call coop_setup_cosmology_from_cosmomc(params, H0/100.d0)\nend subroutine setForH', \
                      r"call setfast(params, CMB)\n call coop_setup_cosmology_from_cosmomc(params)\n call coop_setup_pp()", \
                      r"adjustl(trim(cosmomc_paramnames))" ]


this_index = which_line(param_list, "nnu") 
if(this_index != 0):
    list_param_pattern.append(r'^\s*CMB\%nnu\s*=\s*Params\(\d+\)\s*(\!.*)?$')
    list_param_replace.append(r'CMB%nnu = Params(cosmomc_de_index + cosmomc_de_num_params + ' + str(this_index - index_post_w)  + r' )')



this_index = which_line(param_list, "") 
if(this_index != 0):
    list_param_pattern.append(r'')
    list_param_replace.append(r' = Params(cosmomc_de_index + cosmomc_de_num_params + ' + str(this_index - index_post_w)  + r' )')


this_index = which_line(param_list, "") 
if(this_index != 0):
    list_param_pattern.append(r'')
    list_param_replace.append(r' = Params(cosmomc_de_index + cosmomc_de_num_params + ' + str(this_index - index_post_w)  + r' )')

this_index = which_line(param_list, "yhe") 
if(this_index != 0):
    list_param_pattern.append(r'^\s*CMB\%YHe\s*=\s*params\(\d+\)\s*(\!.*)?$')
    list_param_replace.append(r' CMB%YHe = Params(cosmomc_de_index + cosmomc_de_num_params + ' + str(this_index - index_post_w)  + r' )')


this_index = which_line(param_list, "alpha1") 
if(this_index != 0):
    list_param_pattern.append(r'^\s*CMB\%iso\_cdm\_correlated\s*=\s*params\(\d+\)\s*(\!.*)?$')
    list_param_replace.append(r'CMB%iso_cdm_correlated = Params(cosmomc_de_index + cosmomc_de_num_params + ' + str(this_index - index_post_w)  + r' )')


this_index = which_line(param_list, "deltazrei") 
if(this_index != 0):
    list_param_pattern.append(r'^\s*CMB\%zre\_delta\s*=\s*params\(\d+\)\s*(\!.*)?$')
    list_param_replace.append(r'CMB%zre_delta = Params(cosmomc_de_index + cosmomc_de_num_params + ' + str(this_index - index_post_w)  + r' )')


this_index = which_line(param_list, "Alens") 
if(this_index != 0):
    list_param_pattern.append(r'^\s*CMB\%Alens\s*=\s*params\(\d+\)\s*(\!.*)?$')
    list_param_replace.append(r'CMB%Alens = Params(cosmomc_de_index + cosmomc_de_num_params + ' + str(this_index - index_post_w)  + r' )')


this_index = which_line(param_list, "Alensf") 
if(this_index != 0):
    list_param_pattern.append(r'^\s*CMB\%Alensf\s*=\s*params\(\d+\)\s*(\!.*)?$')
    list_param_replace.append(r'CMB%Alensf = Params(cosmomc_de_index + cosmomc_de_num_params + ' + str(this_index - index_post_w)  + r' )')

this_index = which_line(param_list, "fdm") 
if(this_index != 0):
    list_param_pattern.append( r'^\s*CMB\%fdm\s*=\s*params\(\d+\)\s*(\!.*)?$')
    list_param_replace.append(r'CMB%fdm = Params(cosmomc_de_index + cosmomc_de_num_params + ' + str(this_index - index_post_w)  + r' )')



replace_all("source/CosmologyParameterizations.f90", list_param_pattern, list_param_replace )

replace_first("source/driver.F90", [line_pattern(r'program cosmomc'), line_pattern(r'call ini%open(inputfile)')], [r'program cosmomc \n  use coop_wrapper', r'call Ini%Open(InputFile)\n cosmomc_de_model = Ini%Read_Int("de_model", cosmomc_de_model) \n cosmomc_de_num_params = Ini%Read_Int("de_num_params", cosmomc_de_num_params)\n cosmomc_pp_model = Ini%Read_Int("pp_model", cosmomc_pp_model) \n cosmomc_pp_num_params = Ini%Read_Int("pp_num_params", cosmomc_pp_num_params) \n cosmomc_pp_inflation_consistency = Ini%Read_Logical("inflation_consistency", .true.) \n cosmomc_paramnames = Ini%Read_String("paramnames", .false.) \n if(trim(cosmomc_paramnames).eq."") stop "cannot find the paramnames entry in the ini file" \n cosmomc_de_index = ' + str(index_w) + r'\n cosmomc_de2pp_num_params = ' + str(numhard - 1) + r' - cosmomc_de_index \n cosmomc_pp_num_origin = last_power_index'])

print "************************************"
print "Generating params files"


if(not os.path.exists('paramnames')):
    os.system('mkdir paramnames')
if(iniroot != ""):    
    if(not os.path.exists(inirootname)):    
        os.system('mkdir '+ inirootname )
os.system('cp params_CMB.paramnames paramnames/')
if(num_w_params == 2):
    copy_replace_first("params_CMB.paramnames", "paramnames/params_qcdm.paramnames", [r'^w\s+.+$', r'^wa\s+.+$'], [r"epss          \epsilon_s", r"epsinf           \epsilon_{\infty}\nzetas           \zeta_s\natbyaeq              a_t/a_{\\rm eq}\nQeq                  Q_{\\rm eq} \ndlnQdphi           d\\ln Q/d\\phi" ] )
else:
    copy_replace_first("params_CMB.paramnames", "paramnames/params_qcdm.paramnames", [r'^w\s+.+$'], [r"epss          \epsilon_s\nepsinf           \epsilon_{\infty}\nzetas           \zeta_s\natbyaeq              a_t/a_{\\rm eq}\nQeq                  Q_{\\rm eq} \ndlnQdphi           d\\ln Q/d\\phi" ] )


##lcdm
copy_replace_first(baseini, iniroot+'lcdm.ini', [common_pattern, r'^file_root\s*=.+$', r'^action\s*=.+$',  propose_pattern], [r'\1 \nde_model = 2 \nde_num_params='+str(num_w_params)+r'\npp_model = 0 \npp_num_params = ' + str(index_H0 - index_logA) +r'\nparamnames = paramnames/params_CMB.paramnames', r'file_root = lcdm', r'action = 0', str_propose])

### w0
copy_replace_first_append(baseini, iniroot+'w0.ini', [ r'^file_root\s*=.+$', r'^action\s*=.+$', propose_pattern], [ r'file_root = w0 \nde_model = 2 \nde_num_params='+str(num_w_params)+r'\npp_model = 0 \npp_num_params = ' + str(index_H0 - index_logA) +r'\nparamnames = paramnames/params_CMB.paramnames', r'action = 0', str_propose], r'param[w] = -1 -3 1 0.05 0.05' )

###w0wa
copy_replace_first_append(baseini, iniroot+'w0wa.ini', [ r'^file_root\s*=.+$', r'^action\s*=.+$', propose_pattern], [ r'file_root = w0wa \nde_model = 2 \nde_num_params='+str(num_w_params)+r'\npp_model = 0 \npp_num_params = ' + str(index_H0 - index_logA) +r'\nparamnames = paramnames/params_CMB.paramnames', r'action = 0', str_propose], r'param[w] = -1 -3 1 0.05 0.05 ' + "\n" + r'param[wa] = 0 -3 3 0.1 0.1' )

#qcdm 1 parameter: epss
copy_replace_first_append(baseini, iniroot+'qcdm_1param.ini', [covmat_pattern, r'^file_root\s*=.+$', r'^action\s*=.+$', propose_pattern], [ covmat_repl+r'qcdm_1param.covmat', r'file_root = qcdm_1param \nde_model = 3\nde_num_params=6\npp_model=0 \npp_num_params = ' + str(index_H0 - index_logA) +r'\nparamnames = paramnames/params_qcdm.paramnames', r'action = 0', str_propose], r'param[w] = -1 -1 -1 0 0 '+"\n" + r'param[epss] = 0 -1.5 1.5 0.1 0.1 '+"\n"+r'param[epsinf] = 0 0 0 0 0  '+"\n" + r'param[zetas] = 0 0 0 0 0 '+"\n" +'param[atbyaeq] = 0 0 0 0 0 '+"\n" +'param[Qeq] = 0 0 0 0 0 '+"\n"+'param[dlnQdphi] = 0 0 0 0 0' )


##qcdm 3 parameter: epss, epsinf, zetas
copy_replace_first_append(baseini, iniroot+'qcdm_3param.ini', [covmat_pattern, r'^file_root\s*=.+$', r'^action\s*=.+$',  propose_pattern], [ covmat_repl + r'qcdm_3param.covmat', r'file_root = qcdm_3param \nde_model = 3\nde_num_params=6\npp_model = 0 \npp_num_params = '  + str(index_H0 - index_logA) + r'\nparamnames = paramnames/params_qcdm.paramnames', r'action = 0', str_propose], r'param[w] = -1 -1 -1 0 0 '+"\n"+r'param[epss] = 0 -1.5 1.5 0.1 0.1 '+"\n"+r'param[epsinf] = 0.05 0 1. 0.05 0.05 '+"\n"+r'param[zetas] = 0 -1 1 0.1 0.1 '+"\n"+r'param[atbyaeq] = 0 0 0 0 0 '+"\n"+r'param[Qeq] = 0 0 0 0 0 '+"\n"+r'param[dlnQdphi] = 0 0 0 0 0' )


######################   

 
ppnum = index_H0 - index_logA

for i in range(5, 15):
    ppstr = ''
    for j in range(i):
        ppstr += (r'pp'+ str(j+1) + r'    p_{' + str(j+1) + r'}\n')
    ppstr += r'H0*        H_0'
    copy_replace_first(r'params_CMB.paramnames', r'paramnames/params_scanp' + str(i)+ r'.paramnames',  [ r'^H0\*\s+.+$' ], [ ppstr ])
    copy_replace_first(baseini, iniroot + r'scan_pp'+str(i)+'.ini', [covmat_pattern, common_pattern, r'^file_root\s*=.+$' , r'^action\s*=.+$', r'^compute\_tensors\s*=.+$', r'^(cmb\_dataset\[BICEP2\].*)$', propose_pattern], [covmat_repl + r'scan_pp' + str(i)+r'.covmat', r'DEFAULT(' + batch_dir + '/common_pp.ini) \nde_model = 0 \nde_num_params = 2\npp_model = 1 \npp_num_params = ' + str(ppnum+i) + r'\nparamnames = paramnames/params_scanp' + str(i) + r'.paramnames', r'file_root = scan_pp'+str(i) , r'action = 0', r'compute_tensors = T', r'\#\1',str_propose])
    copy_replace_first_append(baseini, iniroot + r'scan_pp'+str(i)+'_fixrp100.ini', [covmat_pattern, common_pattern, r'^file_root\s*=.+$' , r'^action\s*=.+$', r'^compute\_tensors\s*=.+$', r'^(cmb\_dataset\[BICEP2\].*)$', propose_pattern], [covmat_repl + r'scan_pp' + str(i)+r'_fixrp100.covmat', r'DEFAULT(' + batch_dir + '/common_pp.ini) \nde_model = 0 \nde_num_params = 2\npp_model = 1 \npp_num_params = ' + str(ppnum+i) + r'\nparamnames = paramnames/params_scanp' + str(i) + r'.paramnames', r'file_root = scan_pp'+str(i)+r'_fixrp100' , r'action = 0', r'compute_tensors = T', r'\#\1',str_propose], r'param[r] = 0.1')
    copy_replace_first_append(baseini, iniroot + r'scan_pp'+str(i)+'_fixrp050.ini', [covmat_pattern, common_pattern, r'^file_root\s*=.+$' , r'^action\s*=.+$', r'^compute\_tensors\s*=.+$', r'^(cmb\_dataset\[BICEP2\].*)$', propose_pattern], [covmat_repl + r'scan_pp' + str(i) + r'_fixrp050.covmat', r'DEFAULT(' + batch_dir + '/common_pp.ini) \nde_model = 0 \nde_num_params = 2\npp_model = 1 \npp_num_params = ' + str(ppnum+i) + r'\nparamnames = paramnames/params_scanp' + str(i) + r'.paramnames', r'file_root = scan_pp'+str(i)+r'_fixrp050', r'action = 0', r'compute_tensors = T', r'\#\1',str_propose], r'param[r] = 0.05')
    copy_replace_first_append(baseini, iniroot + r'scan_pp'+str(i)+'_fixrp010.ini', [covmat_pattern, common_pattern, r'^file_root\s*=.+$' , r'^action\s*=.+$', r'^compute\_tensors\s*=.+$', r'^(cmb\_dataset\[BICEP2\].*)$', propose_pattern], [covmat_repl + r'scan_pp'+str(i) + '_fixrp010.covmat', r'DEFAULT(' + batch_dir + '/common_pp.ini) \nde_model = 0 \nde_num_params = 2\npp_model = 1 \npp_num_params = ' + str(ppnum+i) + r'\nparamnames = paramnames/params_scanp' + str(i) + r'.paramnames', r'file_root = scan_pp'+str(i)+r'_fixrp010' , r'action = 0', r'compute_tensors = T', r'\#\1',str_propose], r'param[r] = 0.01')
    copy_replace_first_append(baseini, iniroot + r'scan_pp'+str(i)+'_fixrp001.ini', [covmat_pattern, common_pattern, r'^file_root\s*=.+$' , r'^action\s*=.+$', r'^compute\_tensors\s*=.+$', r'^(cmb\_dataset\[BICEP2\].*)$', propose_pattern], [covmat_repl + r'scan_pp'+str(i) + '_fixrp001.covmat', r'DEFAULT(' + batch_dir + '/common_pp.ini) \nde_model = 0 \nde_num_params = 2\npp_model = 1 \npp_num_params = ' + str(ppnum+i) + r'\nparamnames = paramnames/params_scanp' + str(i) + r'.paramnames', r'file_root = scan_pp'+str(i)+r'_fixrp001' , r'action = 0', r'compute_tensors = T', r'\#\1',str_propose], r'param[r] = 0.001')            
    copy_replace_first(baseini, iniroot + r'scan_pl'+str(i)+'.ini', [covmat_pattern, common_pattern, r'^file_root\s*=.+$' , r'^action\s*=.+$', r'^compute\_tensors\s*=.+$', r'^(cmb\_dataset\[BICEP2\].*)$', propose_pattern], [covmat_repl + r'scan_pl'+str(i) + '.covmat', r'DEFAULT(' + batch_dir + '/common_pp.ini) \nde_model = 0 \nde_num_params = 2\npp_model = 2 \npp_num_params = ' + str(ppnum+i) + r'\nparamnames = paramnames/params_scanp' + str(i) + r'.paramnames', r'file_root = scan_pl'+str(i) , r'action = 0', r'compute_tensors = T', r'\#\1',str_propose])
    copy_replace_first_append(baseini, iniroot + r'scan_pl'+str(i)+'_fixrp100.ini', [covmat_pattern, common_pattern, r'^file_root\s*=.+$' , r'^action\s*=.+$', r'^compute\_tensors\s*=.+$', r'^(cmb\_dataset\[BICEP2\].*)$', propose_pattern], [covmat_repl + r'scan_pl'+str(i) + '_fixrp100.covmat', r'DEFAULT(' + batch_dir + '/common_pp.ini) \nde_model = 0 \nde_num_params = 2\npp_model = 2 \npp_num_params = ' + str(ppnum+i) + r'\nparamnames = paramnames/params_scanp' + str(i) + r'.paramnames', r'file_root = scan_pl'+str(i)+r'fixrp100' , r'action = 0', r'compute_tensors = T', r'\#\1',str_propose], r'param[r] = 0.1')
    

copy_replace_first(common_file, batch_dir + r'/common_pp.ini', [r'^INCLUDE\(params_CMB_defaults\.ini\)\s*$'], [r'INCLUDE(params_CMB_pp.ini)'] )
ppstr = r'param[ns] =  0.96 0.96 0.96 0 0 \n'
for i in range(1, 16):
    ppstr += r'param[pp'+ str(i) + r'] = 0. -1. 1. 0.03 0.03 \n'
copy_replace_first(batch_dir + r'/params_CMB_defaults.ini', batch_dir + r'/params_CMB_pp.ini', [r'^param\[ns\]\s*=.+$', '^param\[r\]\s*=.+$', '^compute_tensors\s*=.+$'], [ ppstr, r'param[r] = 0.12 0. 1. 0.05 0.05', r'compute_tensors = T' ] ) 

copy_replace_first_append(baseini, iniroot+r'scan_bump.ini', [covmat_pattern, r'^file_root\s*=.+$' , r'^action\s*=.+$', r'^compute\_tensors\s*=.+$', r'^(cmb\_dataset\[BICEP2\].*)$', propose_pattern], [covmat_repl + r'scan_bump.covmat', r'file_root = scan_bump \nde_model = 0 \nde_num_params = 2\npp_model = 4 \npp_num_params = ' + str(ppnum+3) + r'\nparamnames = paramnames/params_bump.paramnames' , r'action = 0', r'compute_tensors = T', r'\#\1',str_propose], r'param[bumpamp] = 0. -1. 1. 0.03 0.03 '+"\n"+r'param[bumploc] = -6.3 -6.9 -4.6 0.2 0.2 ' + "\n" + r'param[bumpwidth] = 0.5  0.25 0.75 0.1 0.1' + "\n")

copy_replace_first_append(baseini, iniroot+r'scan_bump_fixrp100.ini', [covmat_pattern, r'^file_root\s*=.+$' , r'^action\s*=.+$', r'^compute\_tensors\s*=.+$', r'^(cmb\_dataset\[BICEP2\].*)$', propose_pattern], [covmat_repl + r'scan_bump_fixrp100.covmat', r'file_root = scan_bump_fixrp100 \nde_model = 0 \nde_num_params = 2\npp_model = 4 \npp_num_params = ' + str(ppnum+3) + r'\nparamnames = paramnames/params_bump.paramnames' , r'action = 0', r'compute_tensors = T', r'\#\1',str_propose], r'param[bumpamp] = 0. -1. 1. 0.03 0.03 '+"\n"+r'param[bumploc] = -6.3 -6.9 -4.6 0.2 0.2 ' + "\n" + r'param[bumpwidth] = 0.5  0.25 0.75 0.1 0.1' + "\n" + r'param[r] = 0.1' + "\n")

copy_replace_first(r'params_CMB.paramnames', r'paramnames/params_bump.paramnames',  [ r'^H0\*\s+.+$' ], [r'bumpamp        A_{\\rm bump} \nbumploc       \\ln k_{\\rm bump} \nbumpwidth     w_{\\rm bump}   \nH0*           H_0'])
