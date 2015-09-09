program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::cosmology
  type(coop_function)::fwp1, fQ
  type(coop_file)::fp
  type(coop_asy)::fig
  COOP_INT::ik
  !!initialize w and alpha functions
  call fwp1%init_polynomial( (/ 0.d0, 0.2d0 /) )
  call fQ%init_polynomial( (/ 0.d0,  0.4d0 /) )
  !!initialize cosmology
  call cosmology%set_coupled_DE_cosmology(Omega_b=0.049d0, Omega_c=0.265d0, h = 0.68d0, tau_re = 0.06d0, As = 2.21d-9, ns = 0.968d0, fwp1 = fwp1, fQ = fQ)
  if(cosmology%h() .eq. 0.d0) stop "Initialization failed; check the input parameters."  
  call cosmology%init_source(0)
  ik = 1
!!$  do while(cosmology%source(0)%k(ik).lt. 0.5d0)
!!$     ik = ik + 1
!!$  enddo
  call cosmology%compute_source_k(cosmology%source(0), ik, do_test_energy_conservation = .true.)
end program test
