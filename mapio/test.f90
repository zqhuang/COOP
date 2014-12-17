program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  implicit none
#include "constants.h"
  integer, parameter::lmax=2000
  COOP_REAL,parameter::fwhm = 15.d0
  type(coop_file)::fp
  integer l, il, i
  COOP_REAL::cls(4, 0:lmax), dtheta, x, als(0:lmax), bls(0:lmax),  ell(0:lmax), sigma, l2cls(4,0:lmax), sigma2, sigma0, sigma1, mu2, mu0, cosbeta

  cls = 0.d0
  sigma = fwhm*coop_sigma_by_fwhm*coop_SI_arcmin
  print*, 1.d0/sigma
  call fp%open("planckbest_lensedtotCls.dat", "r")
  do l=2, lmax
     read(fp%unit, *) il, cls(:, l)
     ell(l)  = l
     l2cls(:,l) = cls(:, l)*coop_2pi*exp(-l*(l+1.d0)*sigma**2)     
     cls(:,l) = l2cls(:,l)/(l*(l+1.d0))
     if(il.ne.l) stop "cl file broken"
  enddo
  call fp%close()
  dtheta = 0.001
  sigma0 = sqrt(sum(Cls(1,0:lmax)*(ell+0.5d0))/coop_2pi)
  sigma1 = sqrt(sum(l2Cls(1,0:lmax)*(ell+0.5d0))/coop_2pi)
  sigma2 = sqrt(sum(l2Cls(1,0:lmax)*(ell+0.5d0)*(ell*(ell+1.d0)))/coop_2pi)
  cosbeta = sigma1**2/sigma0/sigma2
  print*, cosbeta
  call sphere_correlation_init(lmax, als, bls)  
  call fp%open("radial_profile_fwhm15.txt", "w")
  do i = 0, 35
     x = cos(dtheta*i)
     mu2 = sphere_correlation(lmax, l2cls(1, 0:lmax), als, bls, x)/sigma1**2
     write(fp%unit, "(2E16.7)") dtheta*i, sqrt(32.d0/3.d0/coop_pi)*sigma0*mu2*cosbeta
  enddo
  call fp%close()

  
end program test
