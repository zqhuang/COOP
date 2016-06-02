program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::cosmology
  logical::success
  type(coop_file)::fp
  COOP_REAL::kMpc_want, k_want, z_want, tau_want
  type(coop_dictionary)::params
  COOP_INT::ik, output_itau
  type(coop_list_string)::output_variables, copy_var
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
  call coop_dictionary_lookup(params, "variables", output_variables)  
  call coop_dictionary_lookup(params, "kMpc", kMpc_want, 1.1d30)
  if(kMPc_want .gt. 1.d30)then
     call coop_dictionary_lookup(params, "output_redshift", z_want, 1.1d30)
     if(z_want .gt. 1.d30)then
        write(*,*) "Error: for SolvePert you need to set either kMpc or output_redshift in the ini file."
        stop
     endif
  endif

  coop_feedback_level = 3
  call fp%open(output, "w")
  write(fp%unit, "(A)") "#"//trim(params%value("variables"))
  call cosmology%init_from_dictionary(params)
  !!----------------------------------------    
  !!set k/H0  
  call cosmology%init_source(0)
  write(*,*) "H_0 = ", cosmology%h()*100.
  write(*,*) "Omega_m = ", cosmology%omega_m
  if(cosmology%Mpsq0 .ne. 1.d0) write(*,*) "M_*(z=0) = ", sqrt(cosmology%Mpsq0)
  write(*,*)"cosmomc theta (*100) = ",  100.d0 * cosmology%cosmomc_theta()
  write(*,*) "writing the output to file "//trim(output)
  
  if(kMpc_want .gt. 1.d30)then
     tau_want = cosmology%tauofa(1.d0/(1.d0+z_want))
     output_itau =  min(max(1, coop_left_index(cosmology%source(0)%ntau, cosmology%source(0)%tau, tau_want)), cosmology%source(0)%ntau-1)
     if(abs(cosmology%source(0)%tau(output_itau + 1) - tau_want) .lt. abs(cosmology%source(0)%tau(output_itau) - tau_want))then
        output_itau = output_itau + 1
     endif
     write(fp%unit, "(A20, E16.7)") "# fixed  z = ", 1.d0/cosmology%aoftau(cosmology%source(0)%tau(output_itau))-1.d0
     do ik = 1, cosmology%source(0)%nk
        call cosmology%compute_source_k(cosmology%source(0), ik, output = fp%unit, names= output_variables, output_itau = output_itau, success = success)
     enddo
     if(.not. success)stop "Solution blows up exponentially. Model is ruled out."
     
  else
     k_want = kMpc_want/cosmology%H0Mpc()  
     ik = min(max(1, coop_left_index(cosmology%source(0)%nk, cosmology%source(0)%k, k_want)), cosmology%source(0)%nk-1)
     if(abs(cosmology%source(0)%k(ik) - k_want) .gt. abs(cosmology%source(0)%k(ik+1) - k_want))then
        ik = ik + 1
     endif


     write(fp%unit, "(A20, E16.7)") "# fixed  k [Mpc^{-1}] = ", cosmology%source(0)%k(ik)*cosmology%H0Mpc()
     !!--------------solve mode k---------------
     !!output are
     call cosmology%compute_source_k(cosmology%source(0), ik, output = fp%unit, names= output_variables, success = success)
     if(.not. success)write(*,*) "Solution blows up exponentially. Model is ruled out."
  endif
  call fp%close()
end program test
