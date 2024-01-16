program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_REAL:: fnu
  type(coop_arguments)::args
  COOP_REAL,parameter::  sigma0 = 1.d0,sigma2 = 1.d0, xmin = -2.d0, xmax=5.d0, beta_min = 0.2, beta_max = 0.9
  COOP_REAL sigma1
  COOP_INT,parameter::dim = 2, nnu = 80, nbeta = 5, nfit = 4
  COOP_REAL::g(nnu, nbeta), gapp(nnu, nbeta), cosbeta, nu, nuarr(nnu), xarr(nnu), betaarr(nbeta)
  COOP_INT::ic
  COOP_INT::ibeta, inu
  type(coop_file)::fp
  call coop_set_uniform(nnu, xarr, xmin, xmax)
  call coop_set_uniform(nbeta, betaarr, beta_min, beta_max)
  do ibeta = 1, nbeta
     cosbeta = betaarr(ibeta)
     nuarr = xarr/cosbeta
     sigma1 = sqrt(sigma0*sigma2*cosbeta)
     call coop_gaussian_npeak_set_args(args, dim, sigma0, sigma1, sigma2)
     do inu=1, nnu
        nu = nuarr(inu)
        g(inu, ibeta) = log(coop_gaussian_peak_G_slow(nu, args))
        gapp(inu, ibeta) = log(coop_gaussian_peak_G(nu, args))
     enddo
  enddo
  call fp%open("fast.txt", "w")
  do inu = 1,nnu
     write(fp%unit, "("//COOP_STR_OF(nbeta+1)//"E16.7)") xarr(inu), g(inu, :)
  enddo
  call fp%close()
  call fp%open("slow.txt", "w")
  do inu = 1,nnu
     write(fp%unit, "("//COOP_STR_OF(nbeta+1)//"E16.7)")xarr(inu), gapp(inu, :)
  enddo
  call fp%close()

  
end program Test
