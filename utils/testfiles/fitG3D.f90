program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_REAL:: fnu
  type(coop_arguments)::args
  COOP_REAL,parameter::  sigma0 = 1.d0,sigma2 = 1.d0, xmin = -2.d0, xmax=5.d0, beta_min = 0.2, beta_max = 0.95
  COOP_REAL sigma1
  COOP_REAL::bfh_A0 = 1.d0    
  COOP_REAL::bfh_A1 = 0.25d0  
  COOP_REAL::bfh_A2 = 1.d0
  COOP_REAL::bfh_A3 = 1.d0
  COOP_REAL::bfh_B1 = 0.1d0
  COOP_REAL::bfh_B2 = 10.d0
  COOP_REAL::bfh_B3 = 1.d0
  COOP_REAL::bfh_B4 = 1.5d0
  COOP_INT,parameter::dim = 3, nnu = 80, nbeta = 25, nfit = 4
  COOP_REAL::g(nnu, nbeta), gapp(nnu, nbeta), cosbeta, nu, nuarr(nnu), vtest, coef_fit(nfit, nbeta), xarr(nnu), c2d(nfit, nfit), betaarr(nbeta)
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
        g(inu, ibeta) = log(coop_integrate(dgdv, -10.d0, 0.d0, 1.d-9)/args%r(6))+ nu**2/2.d0 - log(coop_2pi)/2.d0
     enddo
     call coop_chebfit(nnu, xarr, g(:, ibeta), nfit, xmin, xmax, coef_fit(:, ibeta))
     print*, coef_fit(:, ibeta)
     do inu=1, nnu
        call coop_chebeval(nfit, xmin, xmax, coef_fit(:, ibeta), xarr(inu), gapp(inu, ibeta))

     enddo
  enddo
  coef_fit= log(abs(coef_fit))
  do ic=1,nfit
     call coop_chebfit(nbeta, betaarr, coef_fit(ic, :),  nfit, beta_min, beta_max, c2d(:,ic))
     print*
     print*, ic
     print*, c2d(:,ic)
     print*
  enddo
  call fp%open("num.txt", "w")
  do inu = 1,nnu
     write(fp%unit, "("//COOP_STR_OF(nbeta+1)//"E16.7)") xarr(inu), g(inu, :)
  enddo
  call fp%close()
  call fp%open("bbks.txt", "w")
  do inu = 1,nnu
     write(fp%unit, "("//COOP_STR_OF(nbeta+1)//"E16.7)")xarr(inu), gapp(inu, :)
  enddo
  call fp%close()

contains

  
  function dGdv(v)
    COOP_REAL dGdv, v
    dGdv = coop_gaussian_npeak_differential(nu, v, args)
  end function dGdv
  

  function gapp_bbks()
    COOP_REAL w, A, B, C1, C2, C3, gapp_bbks
    w = cosbeta*nu
    A = 2.5/(9-5.*cosbeta**2)
    B= 432/sqrt(10.*coop_pi)/(9.-5.*cosbeta*2)**2.5
    C1 = 1.84+1.13*(1-cosbeta**2)**5.72
    C2 = 8.91 + 1.27*exp(6.51*cosbeta**2)
    C3 = 2.58*exp(1.05*cosbeta**2)
    gapp_bbks = (w**3- 3.*w*cosbeta**2+(B*w**2+c1)*exp(-A*w**2))/(1.+c2*exp(-c3*w))
    gapp_bbks = gapp_bbks * (sqrt(coop_2pi)*exp(-nu**2/2.d0))    
  end function gapp_bbks

  function gapp_bfh()
    COOP_REAL::x, gapp_bfh
    x = cosbeta * nu
    gapp_bfh = 0.
  end function gapp_bfh
  
end program Test
