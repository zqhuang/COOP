program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  !!----------------------------------------
  !!output Cls file
  COOP_STRING::output = "cls_coupled_DE_Q2.txt"
  !!----------------------------------------
  COOP_REAL, parameter::hub = 0.676d0 !! h = H_0/100
  COOP_REAL, parameter::ombh2 = 0.022d0 !!Omega_b h^2
  COOP_REAL, parameter::omch2 = 0.12d0   !!Omega_c h^2
  COOP_REAL, parameter::tau_re = 0.08  !!optical depth
  COOP_REAL, parameter::As = 2.219795d-9
  COOP_REAL, parameter::ns = 0.96d0
  
  COOP_REAL, parameter::Omega_b = ombh2/hub**2
  COOP_REAL, parameter::Omega_c = omch2/hub**2

  COOP_REAL, parameter::epsilon_s = 0.25d0
  COOP_REAL, parameter::epsilon_inf = 0.d0
  COOP_REAL, parameter::zeta_s = 0.d0
  COOP_REAL, parameter::beta_s = 6.d0
  !! declare other variables
  type(coop_cosmology_firstorder)::cosmology
  type(coop_function)::fwp1, fQ
  COOP_INT, parameter::lmin = 2, lmax = 2500
  COOP_REAL::Cls(coop_num_Cls, lmin:lmax)
  COOP_REAL::norm, lnorm
  COOP_INT::l
  logical success
  type(coop_file)::fp
#if DO_COUPLED_DE  
  !!----------------------------------------
  !!main code
  !!initialize cosmology
  fwp1 = coop_function_constructor(coop_de_wp1_coupled_quintessence, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args = coop_arguments_constructor( r = (/ 1.d0 - Omega_b - Omega_c , epsilon_s, epsilon_inf, zeta_s , beta_s /) ), name = "DE 1+w")
  
  call fQ%init_polynomial( (/ 0.2d0/) )

  
  call cosmology%set_coupled_DE_cosmology(Omega_b=Omega_b, Omega_c=Omega_c, h = hub, tau_re = tau_re , As = As, ns = ns, fwp1 = fwp1, fQ = fQ)
#else
  Write(*,*) "warning: COUPLED DE disabled; using LCDM model."
  write(*,*) "To enable coupled dark energy, set DARK_ENERGY_MODEL=COUPLED_DE in configure.in and recompile the package."
  call cosmology%set_standard_cosmology(Omega_b=Omega_b, Omega_c=Omega_c, h = hubble, tau_re = tau_re, As = As, ns = ns)  
#endif
  if(cosmology%h() .eq. 0.d0) stop "Initialization failed; check the input parameters."  
  call cosmology%compute_source(0, success = success)
  if(.not. success) stop "Solution blows up exponentially; Model is ruled out."
  call cosmology%source(0)%get_all_cls(lmin, lmax, Cls)
  call fp%open(output,"w")
  write(fp%unit, "(A8, 5A16)") "# ell ", "   TT  ",  "   EE  ",  "   TE   ", "Phi_lens Phi_lens ", " T Phi_lens  "
  norm = cosmology%Tcmb()**2*1.d12
  do l = lmin, lmax
     lnorm = l*(l+1.d0)/coop_2pi*norm
     write(fp%unit, "(I8, 5E16.7)") l, Cls(coop_index_ClTT, l)*lnorm,  Cls(coop_index_ClEE, l)*lnorm,  Cls(coop_index_ClTE, l)*lnorm,  Cls(coop_index_ClLenLen, l)*(l*(l+1.d0))**2*norm, Cls(coop_index_ClTLen, l)*(l*(l+1.d0))**1.5*norm
  enddo
  call fp%close()
  
end program test
