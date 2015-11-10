program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::cosmology
  logical::success
  type(coop_file)::fp
  COOP_REAL::kMpc_want, k_want
  type(coop_dictionary)::params
  COOP_INT::ik
  type(coop_list_string)::output_variables
  type(coop_list_string)::ls
  COOP_STRING::params_file, output

  if(iargc() .lt. 1)then
     write(*,*) "Syntax:"
     write(*,*) "./SolvePert  params.ini"
     stop
  endif
  
  call coop_get_Input(1, params_file)
  call coop_load_dictionary(params_file, params)
  call coop_dictionary_lookup(params, "root", output, default_val="test")
  output = trim(adjustl(output))//"_perturbations.dat"
  call coop_dictionary_lookup(params, "kMpc", kMpc_want)
  call coop_dictionary_lookup(params, "variables", output_variables)



  coop_feedback_level = 3
  call fp%open(output, "w")
  write(fp%unit, "(A)") "#"//trim(params%value("variables"))
  call cosmology%init_from_dictionary(params)
  !!----------------------------------------    
  !!set k/H0  
  call cosmology%init_source(0)
  k_want = kMpc_want/cosmology%H0Mpc()  
  ik = min(max(1, coop_left_index(cosmology%source(0)%nk, cosmology%source(0)%k, k_want)), cosmology%source(0)%nk-1)
  if(abs(cosmology%source(0)%k(ik) - k_want) .gt. abs(cosmology%source(0)%k(ik+1) - k_want))then
     ik = ik + 1
  endif


  write(fp%unit, "(A15, E16.7)") "#  k [Mpc^{-1}] = ", cosmology%source(0)%k(ik)*cosmology%H0Mpc()
  write(*,*) "H_0 = ", cosmology%h()*100. 
  write(*,*)"cosmomc theta (*100) = ",  100.d0 * cosmology%cosmomc_theta()
  write(*,*) "writing the output to file "//trim(output)
  !!--------------solve mode k---------------
  !!output are
  call cosmology%compute_source_k(cosmology%source(0), ik, output = fp%unit, names= output_variables, success = success)
  if(.not. success)write(*,*) "Solution blows up exponentially. Model is ruled out."
  call fp%close()
end program test
