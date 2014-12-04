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
  integer, parameter::lmax=1500
  type(coop_file)::fp
  integer l, il, i
  COOP_REAL::cls(4, 0:lmax), dtheta, x, als(0:lmax), bls(0:lmax), corr, corr0, ell(0:lmax), sigma
  cls = 0.d0
  sigma = 10.d0*coop_sigma_by_fwhm*coop_SI_arcmin
  call fp%open("planckbest_lensedtotCls.dat", "r")
  do l=2, lmax
     read(fp%unit, *) il, cls(:, l)
     ell(l)  = l
     cls(:,l) = cls(:, l)*(coop_2pi/(l*(l+1.d0)))/(1.+l*(l+1.d0)*sigma**2)
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
     write(fp%unit, "(2E16.7)") dtheta*i, corr/corr0*105.24
  enddo
  call fp%close()

  
end program test
