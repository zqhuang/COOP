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
  type(coop_file)::fp
  integer l, il, i
  COOP_REAL::cls(4, 0:lmax), dtheta, x, als(0:lmax), bls(0:lmax), corr, corr0, ell(0:lmax), sigma, l2cls(4,0:lmax), sigma2, sigma0, sigma1
  t
  cls = 0.d0
  sigma = 30.d0*coop_sigma_by_fwhm*coop_SI_arcmin
  call fp%open("planckbest_lensedtotCls.dat", "r")
  do l=2, lmax
     read(fp%unit, *) il, cls(:, l)
     ell(l)  = l
     l2cls(:,l) = cls(:, l)*coop_2pi*exp(-l*(l+1.d0)*sigma**2)     
     cls(:,l) = l2cls(:,l)/(l*(l+1.d0))
     if(il.ne.l) stop "cl file broken"
  enddo
  call fp%close()
  dtheta = 0.7757d-3
  corr0 = sum(Cls(1,0:lmax)*(ell+0.5d0))/coop_2pi
  call sphere_correlation_init(lmax, als, bls)  
  call fp%open("radial_profile.txt", "w")
  do i = 0, 45
     x = cos(dtheta*i)
     corr = sphere_correlation(lmax, cls(1, 0:lmax), als, bls, x)
     write(fp%unit, "(2E16.7)") dtheta*i, corr/corr0*80.
  enddo
  call fp%close()

  
end program test
