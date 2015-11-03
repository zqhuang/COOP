program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  !!This program evolves the perturbations for a given wavenumber k
  !!The default outputs are 7 columns: lna , T00, T00/G00 - 1, T0i, T0i/G0i-1, Phi, Psi
  !!to modify the output, you can change the subroutine "coop_pert_object_print" in pert_object.f90
  
  !!----------------------------------------
  !!wave number k, because COOP uses fixed k arrays, the actual k will be the one that is closest to the following number
  COOP_REAL::kMpc_want = 0.5d0
  COOP_STRING::output = "M0_output.dat"

  
  !!cosmological parameters
  COOP_REAL,parameter::ombM2h2 = 0.02222d0
  COOP_REAL,parameter::omcM2h2 = 0.12d0  !!0.12 LCDM  
  COOP_REAL,parameter::hubble = 0.6783d0  !0.6581d0 !0.7008d0 !0.6801d0 
  COOP_REAL,parameter::tau_re = 0.078d0  !!optical depth2
  COOP_REAL,parameter::As = 2.22d-9   !!amplitude
  COOP_REAL, parameter::ns = 0.9655d0   !!tilt
  COOP_REAL::Omega_b, Omega_c
  !!for EFT Dark Energy I have assumed massless neutrinos, if you want to compare with CAMB/CLASS you need to set mnu = 0

  !!DE background EOS
  COOP_REAL, parameter::w0 = -1.d0
  COOP_REAL, parameter::wa = 0.d0    
  logical::w_is_background = .false.  !!if set to be true, w is defined as the effective background w that gives the same expansion history in GR; otherwise w is defined as p_DE/ rho_DE.
#if DO_EFT_DE  
  !!define the alpha parameters
  COOP_REAL, parameter::alpha_M0 = 0.025d0
  COOP_REAL, parameter::alpha_T0 = 0.d0
  COOP_REAL, parameter::alpha_B0 = 0.01d0
  COOP_REAL, parameter::alpha_K0 = 0.1d0
  COOP_REAL, parameter::alpha_H0 = 0.d0
#elif DO_COUPLED_DE
  COOP_REAL, parameter::Q0 = 0.d0
  COOP_REAL, parameter::Qa = 0.d0  
#endif  
  !!----------------------------------------
  !! declare other variables  
  type(coop_cosmology_firstorder)::cosmology
  type(coop_function)::fwp1, fQ, alphaM, alphaB, alphaK, alphaT, alphaH
  logical::success
  type(coop_file)::fp
  COOP_REAL::M0
  COOP_INT::ik
  !!----------------------------------------
  !!main code

  !!set a higher feedback level
  coop_feedback_level = 3
  call fp%open(output, "w")
  !!initialize w functions as a polynomials of a
  call fwp1%init_polynomial( (/ 1.d0+w0+wa, -wa /) )

#if DO_EFT_DE
  write(fp%unit, "(A)") "#Dark Energy Model = Effective field theory DE"
  alphaM = coop_de_alpha_constructor(alpha_M0, "omega")
  alphaT = coop_de_alpha_constructor(alpha_T0, "omega")
  alphaH = coop_de_alpha_constructor(alpha_H0, "omega")
  alphaB = coop_de_alpha_constructor(alpha_B0, "omega")
  alphaK = coop_de_alpha_constructor(alpha_K0, "omega")
  call coop_EFT_DE_set_Mpsq(alphaM)
#endif
  Omega_b = ombM2h2/coop_Mpsq0/hubble**2
  Omega_c = omcM2h2/coop_Mpsq0/hubble**2

#if DO_EFT_DE  
  !!initialize cosmology
  if(w_is_background)then
     call cosmology%set_EFT_cosmology(Omega_b=omega_b, Omega_c=omega_c, h = hubble, Tcmb = COOP_DEFAULT_TCMB, tau_re = tau_re, As = As, ns = ns, wp1_background = fwp1, alphaM = alphaM, alphaK = alphaK, alphaB= alphaB, alphaH = alphaH, alphaT = alphaT)     
  else
     call cosmology%set_EFT_cosmology(Omega_b=Omega_b, Omega_c=Omega_c, h = hubble, Tcmb = COOP_DEFAULT_TCMB, tau_re = tau_re, As = As, ns = ns, wp1 = fwp1, alphaM = alphaM, alphaK = alphaK, alphaB= alphaB, alphaH = alphaH, alphaT = alphaT)
  endif
#elif DO_COUPLED_DE
  write(fp%unit, "(A)") "#Dark Energy Model = Coupled CDM-DE"
  call fQ%init_polynomial( (/ Q0+Qa, -Qa /) )
  call cosmology%set_coupled_DE_cosmology(Omega_b=Omega_b, Omega_c=Omega_c, h = hubble, tau_re = tau_re, As = As, ns = ns, fwp1 = fwp1, fQ = fQ)
#else
  write(fp%unit, "(A)") "#Dark Energy Model = Lambda (DE w = -1 )"
  call cosmology%set_standard_cosmology(Omega_b=Omega_b, Omega_c=Omega_c, h = hubble, tau_re = tau_re, As = As, ns = ns)
#endif
  !!----------------------------------------    
  !!set k/H0  
  call cosmology%init_source(0)
  ik = 1
  do while(cosmology%source(0)%k(ik).lt. kMpc_want/cosmology%H0Mpc() .and. ik .lt. cosmology%source(0)%nk)
     ik = ik + 1
  enddo
  write(fp%unit, "(A15, E16.7)") "#  k [Mpc^{-1}] = ", cosmology%source(0)%k(ik)*cosmology%H0Mpc()
  write(*,*)"100 theta = ",  100.d0 * cosmology%cosmomc_theta()
  !!--------------solve mode k---------------
  !!output are
  call cosmology%compute_source_k(cosmology%source(0), ik, output = fp%unit, success = success)
  if(.not. success)write(*,*) "Solution blows up exponentially. Model is ruled out."
end program test
