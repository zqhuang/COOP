program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::cosmology
  type(coop_function)::wp1, alphaM, alphaB, alphaK, alphaT, alphaH
  type(coop_file)::fp
  type(coop_asy)::fig
  COOP_INT::ik
  !!initialize w and alpha functions
  call wp1%init_polynomial( (/ 0.d0, 0.15d0 /) )
  call alphaM%init_polynomial( (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /) )
  call alphaK%init_polynomial( (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /) )
  call alphaB%init_polynomial( (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /) )
  call alphaH%init_polynomial( (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /) )
  call alphaT%init_polynomial( (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.1d0 /) )      
  !!initialize cosmology
#if DO_EFT_DE  
  call cosmology%set_EFT_cosmology(Omega_b=0.049d0, Omega_c=0.265d0, h = 0.68d0, tau_re = 0.06d0, As = 2.21d-9, ns = 0.968d0, wp1 = wp1, alphaM = alphaM, alphaK = alphaK, alphaB= alphaB, alphaH = alphaH, alphaT = alphaT)
#else
  call cosmology%set_standard_cosmology(Omega_b=0.049d0, Omega_c=0.265d0, h = 0.68d0, tau_re = 0.06d0, As = 2.21d-9, ns = 0.968d0)
#endif  
  !!set k/H0  
  call cosmology%init_source(0)
  ik = 1
  do while(cosmology%source(0)%k(ik).lt. 2.d0)
     ik = ik + 1
  enddo
  call cosmology%compute_source_k(cosmology%source(0), ik, do_test_energy_conservation = .true.)
end program test
