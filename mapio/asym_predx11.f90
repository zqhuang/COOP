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
  COOP_UNKNOWN_STRING, parameter::color_table = "Planck"
  COOP_REAL, parameter::zmin = 2.
  COOP_REAL, parameter::zmax = 65.
  type(coop_healpix_patch)::patch_s, patch_n, patch_all
  COOP_INT, parameter::mmax = 4
  COOP_INT, parameter::n = 30
  COOP_REAL, parameter::dr = coop_SI_degree*2.d0/n
  type(coop_healpix_maps)::map, mask_s, mask_n, mask_all
  COOP_REAL fr_s(0:n), cov_s(0:n, 0:n), fr_n(0:n), cov_n(0:n, 0:n), diff(0:n), fr_all(0:n), cov_all(0:n, 0:n), err(0:n), chisq
  COOP_INT i, nstep
  type(coop_asy)::fig
  call map%read("predx11/predx11_iqu_smoothed_fwhm15arcmin.fits", nmaps_wanted = 1)
  call mask_all%read("predx11/predx11_imask.fits", nmaps_wanted = 1)
  call mask_s%read("predx11/predx11_imask_south.fits", nmaps_wanted = 1)
  call mask_n%read("predx11/predx11_imask_north.fits", nmaps_wanted = 1)
  call patch_s%init("T", n, dr, mmax = mmax)
  call patch_n%init("T", n, dr, mmax = mmax)
  call patch_all%init("T", n, dr, mmax = mmax)
  call map%stack_with_covariance(patch_s, "spots/predx11_iqu_Tmax_threshold0_fwhm15.txt", coop_healpix_patch_get_fr0, n+1, fr_s, cov_s, mask_s)
  call map%stack_with_covariance(patch_n, "spots/predx11_iqu_Tmax_threshold0_fwhm15.txt", coop_healpix_patch_get_fr0, n+1, fr_n, cov_n, mask_n)
  call map%stack_with_covariance(patch_all, "spots/predx11_iqu_Tmax_threshold0_fwhm15.txt", coop_healpix_patch_get_fr0, n+1, fr_all, cov_all, mask_all)
  call patch_s%plot(1, "predx11_T_on_Tmax_south_fwhm15.txt", caption = "South hemisphere: stacked "//trim(coop_num2str(nint(sum(patch_s%nstack)/patch_s%npix)))//" patches", label = "$T(\mu K)$", headless_vectors = .true., color_table = color_table, zmin = zmin, zmax = zmax)
  call patch_n%plot(1, "predx11_T_on_Tmax_north_fwhm15.txt", caption = "North hemisphere: stacked "//trim(coop_num2str(nint(sum(patch_n%nstack)/patch_n%npix)))//" patches", label = "$T(\mu K)$", headless_vectors = .true., color_table = color_table, zmin = zmin, zmax = zmax)
  call patch_all%plot(1, "predx11_T_on_Tmax_north_fwhm15.txt", caption = "Full sky: stacked "//trim(coop_num2str(nint(sum(patch_n%nstack)/patch_n%npix)))//" patches", label = "$T(\mu K)$", headless_vectors = .true., color_table = color_table, zmin = zmin, zmax = zmax)
  do i=0,n
     err(i) = sqrt(cov_all(i, i)*4./(patch_n%nstack_raw + patch_s%nstack_raw))
     print*, fr_s(i), fr_n(i), fr_all(i),  err(i)
  enddo
  call coop_matsym_inverse(cov_all)
  diff = (fr_s - fr_n)/coop_sqrt2
  chisq = dot_product(diff, matmul(cov_all, diff))*(patch_s%nstack_raw + patch_n%nstack_raw)/2.d0
  print*, "chisq = ", chisq
  print*, "Probability = ", coop_IncompleteGamma((n+1.d0)/2.d0, chisq/2.d0)/Gamma((n+1.d0)/2.d0)


  call fig%open("predx11_radial_profile.txt")
  call fig%init(xlabel = "$r$", ylabel = "$T \,(\mu K)$")
  call  coop_asy_curve(fig, patch_all%r, fr_n, legend="north hemisphere ", color="red", linetype = "solid")
  call coop_asy_curve(fig, patch_all%r, fr_s, legend="south hemisphere", linetype = "dashed")
  nstep = (coop_SI_arcmin*15./patch_all%dr)
  err = err * sqrt(real(nstep))
  do i=nstep/2, n, nstep
     call coop_asy_error_bar(fig, patch_all%r(i), fr_n(i), dy_minus = err(i), dy_plus = err(i), color = "black")
  enddo
  call coop_asy_legend(fig, x = 0.015, y=90.)
  call coop_asy_label(fig, label = "$P(\chi^2 > \chi^2_{\rm measure}) = 5.8\times 10^{-8}$", x = 0.015, y = 40.)
  call fig%close()
  
end program test
