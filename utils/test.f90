program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_REAL:: fnu
  type(coop_arguments)::args
  COOP_REAL,parameter::  sigma0 = 1.d0,sigma2 = 1.d0
  COOP_INT,parameter::dim = 2
  
  COOP_REAL sigma1, cosbeta, nmax, umean, vmean, nu, w0, w2
  cosbeta = 0.8d0
  sigma1 = sqrt(sigma0*sigma2*cosbeta)
  call coop_gaussian_npeak_set_args(args, dim, sigma0, sigma1, sigma2)
  nu  = - 6.d0
  nmax = coop_gaussian_nmax(nu, args)
  umean = coop_gaussian_peak_intu(nu, args)/nmax
  vmean = coop_gaussian_peak_intv(nu, args)/nmax
  print*, "<n_max> = ", nmax, sigma2**2/sigma1**2/8.d0/coop_pi/coop_sqrt3
  print*, "<u> = ", umean, sqrt(32.d0/3/coop_pi)*cosbeta
  print*, "<v> = ", vmean, -sqrt(32.d0/3/coop_pi)
  call coop_gaussian_get_nonoriented_stacking_weights(nu, args, w0, w2)
  print*, "w_0 = ", w0, 0.
  print*, "w_2 = ", w2, -sqrt(32.d0/3/coop_pi)
end program Test
