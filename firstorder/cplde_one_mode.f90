program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  !!cosmological parameters
  COOP_REAL,parameter::ombh2 = 0.022d0
  COOP_REAL,parameter::omch2 = 0.12d0  !!0.12 LCDM  
  COOP_REAL,parameter::hubble = 0.676d0  !!H0/100
  COOP_REAL,parameter::tau_re = 0.08d0  !!optical depth2
  COOP_REAL,parameter::As = 2.22d-9   !!amplitude
  COOP_REAL, parameter::ns = 0.96d0   !!tilt
  COOP_REAL, parameter::Omega_b = ombh2/hubble**2
  COOP_REAL, parameter::Omega_c = omch2/hubble**2
  
  type(coop_cosmology_firstorder)::cosmology
  type(coop_function)::fwp1, fQ
  type(coop_file)::fp
  type(coop_asy)::fig
  COOP_INT::ik
  logical::success
#if DO_COUPLED_DE  
  !!initialize w and alpha functions
  call fwp1%init_polynomial( (/ 0.d0, 0.1d0 /) )
  call fQ%init_polynomial( (/ 0.2d0 /) )
  !!initialize cosmology
  call cosmology%set_coupled_DE_cosmology(Omega_b=Omega_b, Omega_c=Omega_c, h = hubble, tau_re = tau_re, As = As, ns = ns, fwp1 = fwp1, fQ = fQ)
#else
  Write(*,*) "warning: COUPLED DE disabled; using LCDM model."
  write(*,*) "To enable coupled dark energy, set DARK_ENERGY_MODEL=COUPLED_DE in configure.in and recompile the package."
  call cosmology%set_standard_cosmology(Omega_b=Omega_b, Omega_c=Omega_c, h = hubble, tau_re = tau_re, As = As, ns = ns)  
#endif
  if(cosmology%h() .eq. 0.d0) stop "Initialization failed; check the input parameters."  
  call cosmology%init_source(0)
  ik = 1
  do while(cosmology%source(0)%k(ik).lt. 10.d0)
     ik = ik + 1
  enddo
  call cosmology%compute_source_k(cosmology%source(0), ik, do_output = .true., success = success)
  if(.not. success)write(*,*) "Solution blows up exponentially. Model is ruled out."  
end program test
