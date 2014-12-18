program TestNpeak
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_REAL:: fnu
  type(coop_arguments)::args
  COOP_REAL,parameter::  sigma0 = 1.d0,sigma2 = 10000.d0
  COOP_INT,parameter::dim = 2
  
  COOP_REAL sigma1, cosbeta, nmax, umean, vmean, nu, w00, w10, w02, w12, mean(3)
  cosbeta = 0.4d0
  sigma1 = sqrt(sigma0*sigma2*cosbeta)
  call coop_gaussian_npeak_set_args(args, dim, sigma0, sigma1, sigma2)
  nu  = 0.d0
  nmax = coop_gaussian_nmax(nu, args)
  umean = coop_gaussian_peak_intu(nu, args)/nmax
  vmean = coop_gaussian_peak_intv(nu, args)/nmax
  call coop_random_init()
  call coop_gaussian_peak_2Dmax_mean(2, mean, getuv, args)
  print*, "<n_max> = ", nmax
  print*, "<u> = ", umean,  mean(1)
  print*, "<v> = ", vmean, mean(2)
  
  call coop_gaussian_get_nonoriented_stacking_weights(nu, args, w00, w10)
  print*, "w_00 = ", w00
  print*, "w_10 = ", w10
  call coop_gaussian_get_oriented_stacking_weights(nu, args, w00, w10, w02, w12)
  print*, "w_00 = ", w00
  print*, "w_10 = ", w10
  print*, "w_02 = ", w02
  print*, "w_12 = ", w12
  
contains

  subroutine getuv(i,q,u,args,uv,want)
    COOP_REAL::i(0:1),q(0:1),u(0:1)
    type(coop_arguments)::args
    COOP_REAL::uv(3)
    logical want
    want = (i(0).gt.nu)
    if(want)then
       uv(1:2) = i
       uv(3) = sqrt(q(0)**2+q(1)**2)
    endif
  end subroutine getuv
  
end program TestNpeak
