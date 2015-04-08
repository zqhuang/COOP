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


def function_pattern(func, var=''):
    if var == '':
        return r'function\s+' + func + r'(\s*\(\s*\))?\s*(\!.*)?\n(.*\n)+\s*end\s+function\s+' + func        
    else:
        return r'function\s+' + func + r'\s*\(\s*' + re.sub(r'\,', r'\s*\,\s*', var) + r'\s*\)\s*(\!.*)?\n(.*\n)+\s*end\s+function\s+' + func

def subroutine_pattern(func, var=''):
    if var == '':
        return r'subroutine\s+' + func + r'(\s*\(\s*\))?\s*(\!.*)?\n(.*\n)+\s*end\s+subroutine\s+' + func        
    else:
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


replace_first("camb/camb.f90", [r'^\s*subroutine\s*CAMB_GetResults\s*\(.*\)\s*(!.*)?$', r'if\s*\(\s*Params%WantCls\s*\.and\.\s*Params%WantScalars\s*\)\s*then'], [r'    subroutine CAMB_GetResults(Params, error)\n   use coop_wrapper', r'call coop_global_cosmology_prepare_from_camb() \n     if (Params%WantCls .and. Params%WantScalars) then'])

replace_all("camb/equations_ppf.f90", \
            [ line_pattern(r"module lambdageneral"),  \
              r'is_cosmological_constant\s*=.+', \
              function_pattern('w_de','a'), \
              function_pattern('grho_de', 'a')],
            [ r'#include "constants.h"\nmodule LambdaGeneral \nuse coop_wrapper', \
              r'is_cosmological_constant = .false.', \
              r'function w_de(a)\n real(dl) w_de, a\n w_de = coop_global_de%weffofa(COOP_REAL_OF(a)) \n end function w_de', \
              r'function grho_de(a)\n real(dl) grho_de, a\n   grho_de = grhov * coop_global_cosmology_DE_rhoa4_ratio_eff(a) \n end function grho_de'] )


replace_first("camb/cmbmain.f90", [subroutine_pattern('DoSourcek', 'Ev,q_ix'), r'call\s+GetTransfer\s*\(\s*Ev\s*\,\s*tau\)'], [r'subroutine DoSourcek(Ev, q_ix) \n use coop_wrapper \n#include "constants.h"\n integer q_ix \n   type(EvolutionVars) EV  \n if(CP%WantScalars .and. global_error_flag == 0) call COOP_COSMO%camb_dosourcek(0, q_ix, Evolve_q%points(q_ix), TimeSteps%points(2:TimeSteps%npoints), src(:,:,2:TimeSteps%npoints), CP%transfer%num_redshifts, tautf(1:CP%transfer%num_redshifts), MT%TransferData) \n if(CP%WantTensors .and. global_error_flag == 0) call COOP_COSMO%camb_dosourcek(2, q_ix, Evolve_q%points(q_ix), TimeSteps%points(2:TimeSteps%npoints), src(:,:,2:timesteps%npoints)) \n end subroutine DoSourcek', r'call COOP_COSMO%camb_gettransfer(q_ix, EV%q, CP%transfer%num_redshifts, tautf(1:CP%transfer%num_redshifts), MT%TransferData)'] )

#replace_all("source/wl.f90", [r'useweyl\s*\=.*'], [ r'useweyl = .true.'])

#replace_all("source/CosmologyTypes.f90", [r'use\_Weylpower\s*\=.*'], [ r'use_weylpower = .true.'])

replace_first("source/settings.f90", [line_pattern(r'integer,parameter::max_theory_params=\d+')], [r'integer, parameter:: max_theory_params = 50 \n character(LEN=256)::cosmomc_paramnames = "'+default_params+'" '])

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
                      r"(\'|\")" + re.escape(default_params) + r"(\'|\")", \
                      r'^\s*(real.*\:\:\s*use_min_zre\s*=.*)$', \
                      r'^\s*if\s*\(\s*CMB\%zre\s*\<\s*this\%Use\_min\_zre\)\s*return\s*(\!.*)?$' ]



list_param_replace = [r'#include "constants.h"\nmodule CosmologyParameterizations\n use coop_wrapper', \
                      r'H0_min = 50.', \
                      r'H0_max = 90.', \
                      r'this%H0_min = 50.', \
                      r'this%H0_max = 90. \n call Ini%read("use_max_zre", this%use_max_zre)', \
                      r'call this%SetTheoryParameterNumbers(cosmomc_de_index + cosmomc_de_num_params+cosmomc_de2pp_num_params-1, cosmomc_pp_num_params)', \
                      r'call coop_setup_cosmology_from_cosmomc(params, H0/100.d0)\nend subroutine setForH', \
                      r"         call setfast(params, CMB)\n    call coop_setup_cosmology_from_cosmomc(params)\n " , \
                      r"adjustl(trim(cosmomc_paramnames))", \
                      r'\1\n real(mcp)::use_max_zre = 20.', \
                      r' if (CMB%zre < this%Use_min_zre .or. CMB%zre > this%use_max_zre) return \n ']


this_index = which_line(param_list, "nnu") 
if(this_index != 0):
    list_param_pattern.append(r'^\s*CMB\%nnu\s*=\s*Params\(\d+\)\s*(\!.*)?$')
    list_param_replace.append(r'CMB%nnu = Params(cosmomc_de_index + cosmomc_de_num_params + ' + str(this_index - index_post_w)  + r' )')


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

if(num_w_params == 2):
    copy_replace_first(default_params, "paramnames/params_coupledDE.paramnames", [r'^w\s+.+$', r'^wa\s+.+$'], [r'Qde          Q_{\\rm DE}', r'nde           n_{\rm DE}\ndlnQdphi             d\\ln Q/d\\phi           \ndUdphi                  dU/d\\phi \nd2Udphi2           d^2U/d\\phi^2' ] )
else:
    copy_replace_first(default_params, "paramnames/params_coupledDE.paramnames", [r'^w\s+.+$'], [r'Qde          Q_{\\rm DE} \nnde           n_{\rm DE}\ndlnQdphi             d\\ln Q/d\\phi\ndUdphi                  dU/d\\phi\nd2Udphi2           d^2U/d\\phi^2'])




##coupled DE 5 parameters: Q, n, dQ/dphi,  dU/dphi, d^2U/dphi^2
mycovmat = ''   # 'coupledDE.covmat'
copy_replace_first_append(baseini, iniroot+'coupledDE.ini', [covmat_pattern, r'^file_root\s*=.+$', r'^action\s*=.+$',  propose_pattern], [ covmat_repl + mycovmat, r'file_root = coupledDE \nde_model = 4\nde_num_params=5\npp_model = 0 \npp_num_params = '  + str(index_H0 - index_logA) + r'\nparamnames = paramnames/params_coupledDE.paramnames', r'action = 4', str_propose], r'param[Qde] = 0. 0. 0.5 0.05 0.05 '+"\n"+r'param[nde] = 0'+"\n"+r'param[dlnQdphi] = 0'+"\n"+r'param[dUdphi] = 0 '+"\n"+r'param[d2Udphi2] = 0' +"\n" + r'use_min_zre = 6.' + "\n" + r'use_max_zre = 15.' + "\n\n" )


