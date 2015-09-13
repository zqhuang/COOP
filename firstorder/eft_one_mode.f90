program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  !!----------------------------------------
  !!wave number k, because COOP uses fixed k arrays, the actual k will be the one that is closest to the following number
  COOP_REAL::kMpc_want = 0.2d0

  !!background EOS
  COOP_REAL, parameter::w0 = -1.d0
  COOP_REAL, parameter::wa = 0.d0    
  !!define the alpha parameters
  COOP_REAL, parameter::alpha_M0 = 0.d0
  COOP_REAL, parameter::alpha_T0 = 0.d0
  COOP_REAL, parameter::alpha_B0 = 0.d0
  COOP_REAL, parameter::alpha_K0 = 0.d0
  COOP_REAL, parameter::alpha_H0 = 0.d0
  !!----------------------------------------
  !! declare other variables  
  type(coop_cosmology_firstorder)::cosmology
  type(coop_function)::wp1, alphaM, alphaB, alphaK, alphaT, alphaH
  logical::success
  COOP_INT::ik
  !!----------------------------------------
  !!main code
#if DO_EFT_DE
  
  !!----------------------------------------  
  !!initialize w functions as a polynomials of a
  call wp1%init_polynomial( (/ 1.d0+w0+wa, -wa /) )

  
  !!----------------------------------------  
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

  
  !!----------------------------------------    
  !!initialize cosmology
  call cosmology%set_EFT_cosmology(Omega_b=0.049d0, Omega_c=0.265d0, h = 0.68d0, tau_re = 0.06d0, As = 2.21d-9, ns = 0.968d0, wp1 = wp1, alphaM = alphaM, alphaK = alphaK, alphaB= alphaB, alphaH = alphaH, alphaT = alphaT)
#else
  Write(*,*) "warning: EFT DE disabled; using LCDM model."
  write(*,*) "To enable EFT dark energy, set DARK_ENERGY_MODEL=EFT in configure.in and recompile the package."
  call cosmology%set_standard_cosmology(Omega_b=0.049d0, Omega_c=0.265d0, h = 0.68d0, tau_re = 0.06d0, As = 2.21d-9, ns = 0.968d0)
#endif
  
  !!----------------------------------------    
  !!set k/H0  
  call cosmology%init_source(0)
  ik = 1
  do while(cosmology%source(0)%k(ik).lt. kMpc_want/cosmology%H0Mpc() .and. ik .lt. cosmology%source(0)%nk)
     ik = ik + 1
  enddo
  write(*,*) "k [Mpc^{-1}] = ", cosmology%source(0)%k(ik)*cosmology%H0Mpc()
  !!--------------solve mode k---------------
  !!output are
  write(*,"(5A16)") "# ln a",    " T00 ",  " T00/G00 - 1 ",   "T0i, T0i/G0i - 1"
  
  call cosmology%compute_source_k(cosmology%source(0), ik, do_test_energy_conservation = .true., success = success)
  if(.not. success)write(*,*) "Solution blows up exponentially. Model is ruled out."
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
