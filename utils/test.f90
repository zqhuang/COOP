program TestNpeak
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_REAL:: fnu
  type(coop_arguments)::args
  COOP_REAL,parameter::  sigma0 = 1.d0,sigma2 = 5.d0
  COOP_INT,parameter::dim = 2
  
  COOP_REAL sigma1, cosbeta, nmax, umean, vmean, nu,  mean(2), weights(4), pmean(2), twonu2
  cosbeta = 0.6d0
  sigma1 = sqrt(sigma0*sigma2*cosbeta)
  call coop_gaussian_npeak_set_args(args, dim, sigma0, sigma1, sigma2)
  nu  = 1.d0
  
  twonu2  = 2.d0*nu**2
  
  nmax = coop_gaussian_nmax(nu, args)
  umean = coop_gaussian_peak_intu(nu, args)/nmax
  vmean = coop_gaussian_peak_intv(nu, args)/nmax
  call coop_random_init()
  call coop_gaussian_peak_2Dmax_mean(2, mean, getuv, args)
  print*, "<n_max> = ", nmax
  print*, "<u> = ", umean,  mean(1)
  print*, "<v> = ", vmean, mean(2)
  call coop_gaussian_peak_Pmax_mean(2, pmean, getp, args)
  print*, "<Pmax>, <Laplacian Pmax> = ", pmean
  call coop_gaussian_get_nonoriented_stacking_weights(nu, args, weights)
  write(*,"(4G14.5)") weights
  call coop_gaussian_get_oriented_stacking_weights(nu, args, weights)
  write(*,"(4G14.5)") weights

  
contains

  subroutine getp(q, u, qx, qy, qq, uu, qu, uq, args, s, want)
    COOP_REAL:: q(0:1), u(0:1), qx, qy, qq, uu, qu, uq
    type(coop_arguments)::args
    COOP_REAL::s(2)
    logical want
    want =  (q(0)**2+u(0)**2.gt. twonu2)
    if(want)then
       s(1) = sqrt(q(0)**2+u(0)**2)
       s(2) = (q(1)*q(0)+u(1)*u(0))/s(1)
    endif
  end subroutine getp

  subroutine getuv(i,q,u,args,uv,want)
    COOP_REAL::i(0:1),q(0:1),u(0:1)
    type(coop_arguments)::args
    COOP_REAL::uv(2)
    logical want
    want = (i(0).gt.nu)
    if(want)then
       uv(1:2) = i
    endif
  end subroutine getuv
  
end program TestNpeak
