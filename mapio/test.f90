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
  COOP_REAL,parameter::nu = 1.d0
  type(coop_file)::fp0, fp2, fp
  integer l, il, i
  type(coop_arguments)::args
  
  COOP_REAL::cls(4, 0:lmax), dtheta, x, als(0:lmax), bls(0:lmax),  ell(0:lmax), sigma, l2cls(4,0:lmax), sigma2, sigma0, sigma1, cosbeta, w00, w10, w02,w12, c00, c10, c02, c12, j2, theta, omega
  call coop_random_init()
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
  call coop_gaussian_get_oriented_stacking_weights(nu, args, w00,w10,w02,w12)
  print*, w00*sigma0, w10*sigma2,w02*sigma0,w12*sigma2
  call sphere_correlation_init(lmax, als, bls)  
  call fp0%open("T0_fwhm"//COOP_STR_OF(nint(fwhm))//".txt", "w")
  call fp2%open("T2_fwhm"//COOP_STR_OF(nint(fwhm))//".txt", "w")  
  do i = 0, n
     theta = dtheta*i
     omega = 2.d0*sin(theta/2.d0)
     x = cos(theta)
     c00 = sphere_correlation(lmax, cls(1, 0:lmax), als, bls, x)
     c10 = -sphere_correlation(lmax, l2cls(1, 0:lmax), als, bls, x)
     c02 = 0.d0
     c12 = 0.d0     
     do l = 0, lmax
        j2 = coop_bessj(2, (l+0.5d0)*theta)
        c02 = c02 + (l+0.5d0)*j2*cls(1,l)
        c12 = c12 - (l+0.5d0)*j2*l2cls(1,l)        
     enddo
     c02 = c02/coop_2pi
     c12 = c12/coop_2pi
     write(fp0%unit, "(4E16.7)") omega, w00*c00+w10*c10
     write(fp2%unit, "(4E16.7)") omega, w02*c02+w12*c12     
  enddo
  call fp0%close()
  call fp2%close()

  
end program test
