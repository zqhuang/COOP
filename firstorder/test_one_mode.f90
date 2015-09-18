program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  !!----------------------------------------
  !!wave number k, because COOP uses fixed k arrays, the actual k will be the one that is closest to the following number
  COOP_REAL::kMpc_want = 0.2d0

  !!----------------------------------------
  !! declare other variables  
  type(coop_cosmology_firstorder)::cosmology
  type(coop_function)::wp1, alphaM, alphaB, alphaK, alphaT, alphaH
  logical::success
  COOP_INT::ik
  !!----------------------------------------
  !!main code
  !!----------------------------------------    
  !!initialize cosmology
  call cosmology%set_standard_cosmology(Omega_b=0.049d0, Omega_c=0.265d0, h = 0.68d0, tau_re = 0.06d0, As = 2.21d-9, ns = 0.968d0)
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
  
  call cosmology%compute_source_k(cosmology%source(0), ik, do_output = .true., success = success)
  if(.not. success)write(*,*) "Solution blows up exponentially. Model is ruled out."
  
end program test
