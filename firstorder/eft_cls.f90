program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  !!----------------------------------------
  !!output Cls file
  COOP_STRING::output = "cls_scalar_output.txt"
  !!background EOS
  COOP_REAL, parameter::w0 = -1.d0
  COOP_REAL, parameter::wa = 0.d0    
  !!define the alpha parameters  
  COOP_REAL, parameter::alpha_M0 = 0.d0
  COOP_REAL, parameter::alpha_T0 = 0.d0
  COOP_REAL, parameter::alpha_K0 = 0.d0
  COOP_REAL, parameter::alpha_B0 = 0.d0
  COOP_REAL, parameter::alpha_H0 = 0.d0
  !!----------------------------------------
  COOP_REAL, parameter::hub = 0.676d0 !! h = H_0/100
  COOP_REAL, parameter::ombh2 = 0.022d0 !!Omega_b h^2
  COOP_REAL, parameter::omch2 = 0.12d0   !!Omega_c h^2
  COOP_REAL, parameter::tau_re = 0.08  !!optical depth
  COOP_REAL, parameter::As = 2.219795d-9
  COOP_REAL, parameter::ns = 0.96d0
  !! declare other variables
  type(coop_cosmology_firstorder)::cosmology
  type(coop_function)::wp1, alphaM, alphaB, alphaK, alphaT, alphaH
  COOP_INT, parameter::lmin = 2, lmax = 2850
  COOP_REAL::Cls(coop_num_Cls, lmin:lmax)
  COOP_REAL::norm, lnorm
  COOP_INT::l
  logical success
  type(coop_file)::fp
  !!----------------------------------------
  !!main code
#if DO_EFT_DE
  
  !!initialize w as a polynomial of a
  call wp1%init_polynomial( (/ 1+w0+wa, -wa /) )
  
  !!initialize alpha functions as  alpha_X(a) = alpha_X0 H_0^2/H(a)^2, where H(a) is LCDM Hubble 
  call generate_function(alpha_M0, alphaM)
  call generate_function(alpha_T0, alphaT)
  call generate_function(alpha_H0, alphaH)
  call generate_function(alpha_B0, alphaB)
  call generate_function(alpha_K0, alphaK)
  !!------- uncomment the following lines if you want to use polynomials, too
  !!call alphaM%init_polynomial( (/ 0.d0, 0.d0, 0.1d0 /) )
  !!call alphaK%init_polynomial( (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /) )
  !!call alphaB%init_polynomial( (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /) )
  !!call alphaH%init_polynomial( (/ 0.d0, 0.d0,  0.d0, 0.d0 /) )
  !!call alphaT%init_polynomial( (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /) )

  !!initialize cosmology
  call cosmology%set_EFT_cosmology(Omega_b=ombh2/hub**2, Omega_c=omch2/hub**2, h = hub, tau_re = tau_re, As = As, ns = ns, wp1 = wp1, alphaM = alphaM, alphaK = alphaK, alphaB= alphaB, alphaH = alphaH, alphaT = alphaT)
  call cosmology%compute_source(0, success = success)
  if(.not. success) stop "Solution blows up exponentially; Model is ruled out."
  call cosmology%source(0)%get_all_cls(lmin, lmax, Cls)
  call fp%open(output,"w")
!  write(fp%unit, "(A8, 5A16)") "# ell ", "   TT  ",  "   EE  ",  "   TE   ", "Phi_lens Phi_lens ", " T Phi_lens  "
  do l = lmin, lmax-100
     norm =cosmology%Tcmb()**2*1.d12
     lnorm =  l*(l+1.d0)/coop_2pi*norm
     write(fp%unit, "(I8, 5E16.7)") l, Cls(coop_index_ClTT, l)*lnorm,  Cls(coop_index_ClEE, l)*lnorm,  Cls(coop_index_ClTE, l)*lnorm,  Cls(coop_index_ClLenLen, l)*(l*(l+1.d0))**2*norm, Cls(coop_index_ClTLen, l)*(l*(l+1.d0))**1.5*norm
  enddo
  call fp%close()
#else
  write(*,*) "EFT dark energy disabled."
  write(*,*) "To enable EFT dark energy, set DARK_ENERGY_MODEL=EFT in configure.in and recompile the package."
#endif

contains

  subroutine generate_function(alpha0, f)
    COOP_REAL::alpha0
    type(coop_function)::f    
    COOP_REAL,parameter::omega_m0= 0.3
    COOP_REAL,parameter::omega_r0 = 8.d-5
    COOP_INT, parameter::n = 8192
    COOP_REAL::a(n), alpha(n)
    COOP_INT::i
    call coop_set_uniform(n, a, log(coop_min_scale_factor), log(coop_scale_factor_today))
    a = exp(a)
    alpha = alpha0/(1.d0-omega_m0-omega_r0 + omega_m0/a**3 + omega_r0/a**4)
    call f%init(n = n, xmin = a(1), xmax = a(n), f = alpha, xlog = .true., ylog = .false.)
  end subroutine generate_function
  
end program test
