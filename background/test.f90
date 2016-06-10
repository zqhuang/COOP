program test
  use coop_background_mod
  use coop_wrapper_utils
#include "constants.h"
  type(coop_function)::fofR, Vofphi
  COOP_REAL,parameter::Omega_m = 0.3
  COOP_REAL,parameter::Omega_Lambda = 1.d0 - Omega_m
  COOP_REAL::Lambda,  n, c2, phi
  COOP_INT::i
  !! action = 1/2 \int [R - 2 Lambda + f(R) ]...
  !!f(R) =   2 Lambda /(c2 * R^n + 1)
  Lambda = 3.d0*Omega_Lambda
  n = 4.d0
  c2 = 12.*Omega_Lambda/Omega_m ! 1.d6*n/Lambda**n

  call fofR%init_rational( c_up =(/ 2*Lambda /), alpha_up = (/ 0.d0 /), c_down = (/ c2, 1.d0 /),  alpha_down = (/ n, 0.d0 /) )
  call coop_convert_fofR_to_Vofphi(fofR, 2.1d0, Vofphi)
  do i = 1, Vofphi%n
     phi = exp(Vofphi%f1(i))
     write(*,*) phi, Vofphi%eval(phi), Vofphi%derivative(phi)
  enddo
end program test
