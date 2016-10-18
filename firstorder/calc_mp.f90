program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  COOP_STRING::params_file, output
  type(coop_dictionary)::params  
  COOP_REAL::redshift
  COOP_INT, parameter::nk = 256
  COOP_INT::ik
  COOP_REAL k(nk), matterPk(nk), khMpc(nk), PsiPk(nk), PhiPk(nk), khMpc_min, khMpc_max
  type(coop_cosmology_firstorder)::cosmology
  logical::success
  type(coop_file)::fp
  !!----------------------------------------
  !!----------------------------------------
  if(iargc() .lt. 1)then
     write(*,*) "Syntax:"
     write(*,*) "./CalcMP  params.ini"
     stop
  endif
  coop_feedback_level = 4  
  coop_source_k_index = 0.25 !!sample more low k
  call coop_get_Input(1, params_file)
  call coop_load_dictionary(params_file, params)
  call cosmology%init_from_dictionary(params)
  call coop_dictionary_lookup(params, "root", output, default_val="test")
  call coop_dictionary_lookup(params, "output_redshift", redshift, default_val = 0.d0)
  call coop_dictionary_lookup(params, "khMpc_min", khMpc_min, default_val = 1.d-4)
  call coop_dictionary_lookup(params, "khMpc_max", khMpc_max, default_val = 0.5d0)
  write(*,"(A, F12.3)") "Computing power spectrum at redshift z = ", redshift
  call cosmology%init_source(0)
  call cosmology%compute_source(0, success = success)  
  write(*,*) "Sigma_8 = ", cosmology%sigma_8
  if(.not. success) stop "Unhealthy model: perturbations blow up exponentially."
  call coop_set_uniform(nk, khMpc, khMpc_min, khMpc_max, logscale = .true.)
  k = khMpc *cosmology%h()/cosmology%H0Mpc()  !!k

  !!compute k^3 |\delta_k|^2 /(2pi^2)  
  call cosmology%get_Matter_power(z=redshift, nk = nk, k = k, Pk = matterPk)
  !!compute k^3 [ k^2/a^2 Phi_k /(3/2H^2\Omega_m) ]^2/(2pi^2)
  call cosmology%get_Phi_power(z=redshift, nk = nk, k = k, Pk = PhiPk)
  !!compute k^3 [ k^2/a^2 Psi_k /(3/2H^2\Omega_m) ]^2/(2pi^2)  
  call cosmology%get_Psi_power(z=redshift, nk = nk, k = k, Pk = PsiPk)
  !!for GR matterPK = PsiPK = PhiPK
  

  matterPk = matterPk * (2.d0*coop_pi**2)/khMpc**3  !!obtain |\delta_k|^2 in unit of (Mpc/h)^3
  PsiPk = PsiPk * (2.d0*coop_pi**2)/khMpc**3  !!obtain |\delta_k|^2 in unit of (Mpc/h)^3
  PhiPk = PhiPk * (2.d0*coop_pi**2)/khMpc**3  !!obtain |\delta_k|^2 in unit of (Mpc/h)^3    
  write(*,*) "Saving the matter power spectrum in "//trim(output)//"_matterpower.dat"
  call fp%open(trim(output)//"_matterpower.dat","w")
  write(fp%unit, "(4A22)") "#k [h/Mpc]   ",  " CDM power[(Mpc/h)^3]", " Phi power[(Mpc/h)^3]", " Psi power[(Mpc/h)^3]"
  do ik=1, nk
     write(fp%unit, "(4E22.7)") khMpc(ik), matterPk(ik), PhiPk(ik), PsiPk(ik)
  enddo
  call fp%close()
  
end program test
