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
  integer, parameter::lmax=2000, n=50
  COOP_REAL,parameter::fwhm = 15.d0
  COOP_REAL,parameter::nu = 0
  type(coop_file)::fp
  integer l, il, i
  type(coop_arguments)::args
  
  COOP_REAL::cls(4, 0:lmax), dtheta, x, als(0:lmax), bls(0:lmax),  ell(0:lmax), sigma, l2cls(4,0:lmax), sigma2, sigma0, sigma1, mu2, mu0, cosbeta, w0, w2

  cls = 0.d0
  sigma = fwhm*coop_sigma_by_fwhm*coop_SI_arcmin
  call fp%open("planckbest_lensedtotCls.dat", "r")
  do l=2, lmax
     read(fp%unit, *) il, cls(:, l)
     ell(l)  = l
     l2cls(:,l) = cls(:, l)*coop_2pi*exp(-l*(l+1.d0)*sigma**2)     
     cls(:,l) = l2cls(:,l)/(l*(l+1.d0))
     if(il.ne.l) stop "cl file broken"
  enddo
  call fp%close()
  dtheta = 2.d0*coop_SI_degree/n
  sigma0 = sqrt(sum(Cls(1,0:lmax)*(ell+0.5d0))/coop_2pi)
  sigma1 = sqrt(sum(l2Cls(1,0:lmax)*(ell+0.5d0))/coop_2pi)
  sigma2 = sqrt(sum(l2Cls(1,0:lmax)*(ell+0.5d0)*(ell*(ell+1.d0)))/coop_2pi)
  cosbeta = sigma1**2/sigma0/sigma2
  print*, "gamma = ", cosbeta
  print*, "sigma_0 = ", sigma0
  print*, "sigma_2 = ", sigma2  
  call coop_gaussian_npeak_set_args(args, 2, sigma0, sigma1, sigma2)  
  call coop_gaussian_get_nonoriented_stacking_weights(nu, args, w0, w2)
  print*, w0, w2
  call sphere_correlation_init(lmax, als, bls)  
  call fp%open("radial_profile_fwhm"//COOP_STR_OF(nint(fwhm))//".txt", "w")
  do i = 0, n
     x = cos(dtheta*i)
     mu2 = -sphere_correlation(lmax, l2cls(1, 0:lmax), als, bls, x)
     mu0 = sphere_correlation(lmax, cls(1, 0:lmax), als, bls, x)
     write(fp%unit, "(4E16.7)") dtheta*i, mu0*w0 + mu2*w2
  enddo
  call fp%close()

  
end program test
