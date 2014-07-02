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
  type(coop_healpix_patch)::patch
  COOP_INT, parameter::mmax = 4
  COOP_INT, parameter::n = 30
  COOP_REAL, parameter::dr = coop_SI_degree*2.d0/n
  type(coop_healpix_maps)::map, mask
  COOP_REAL fr(0:n), cov(0:n, 0:n)
  COOP_INT i

  call map%read("predx11/predx11_iqu_smoothed_fwhm15arcmin.fits", nmaps_wanted = 1)
!  call map%smooth(15.d0*coop_SI_arcmin)
  call mask%read("predx11/predx11_imask.fits", nmaps_wanted = 1)
  call patch%init("T", n, dr, mmax = mmax)
  call map%stack_with_covariance(patch, "spots/predx11_iqu_Tmax_threshold0_fwhm15.txt", coop_healpix_patch_get_fr0, n+1, fr, cov, mask)
!  call map%stack(patch, "spots/predx11_iqu_Tmax_threshold0_fwhm15.txt",  mask)
  write(*,*) "stacked "//trim(coop_num2str(nint(sum(patch%nstack)/patch%npix)))//" patches"

  call patch%plot(1, "test.txt", caption = "stacked "//trim(coop_num2str(nint(sum(patch%nstack)/patch%npix)))//" patches", label = "$T(\mu K)$", headless_vectors = .true.)
  do i=0,n
     print*, fr(i), sqrt(cov(i,i))
  enddo
  
end program test
