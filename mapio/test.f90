program shells
  use coop_wrapper_firstorder
  use coop_zeta3d_mod
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::cosmology
  type(coop_function)::deltaN
  COOP_REAL,parameter::sigma_zeta = 4.6d-5
  COOP_INT,parameter::lmax = 250
  COOP_INT,parameter::nside = 256
  COOP_REAL,parameter::fwhm_arcmin = 120.d0
  COOP_REAL::sigma_chi = 1.d-7
  COOP_REAL::mean_chi = 5.8d-6
  COOP_REAL::gp_width = 3.
  COOP_REAL::gp_A = 2.e-4
  COOP_REAL::gp_mean = 12.
  COOP_SINGLE::resolution = 1. !!In unit of h^{-1}Mpc
  COOP_INT:: i, j

  call coop_file_load_function("deltaN-LUT-1.875", 1, 2, deltaN, .false.)
  call cosmology%Set_Planck_bestfit()
  call cosmology%compute_source(0)
  call coop_feedback("source done")
#ifdef BILLIARDS
  write(*,*) "<chi> = ", mean_chi
  call coop_zeta3d_generate_cmb( cosmology, fnl, lmax, nside, "zetaproj/testmap", "zetaproj/gp_meanchi"//COOP_STR_OF(nint(mean_chi*1.e7))//"_sigmachi"//COOP_STR_OF(nint(sigma_chi*1.e7)) , fwhm_arcmin = fwhm_arcmin)
#else
  call coop_zeta3d_generate_cmb( cosmology, fnl, lmax, nside, "zetaproj/testmap", "zetaproj/gs_amp"//COOP_STR_OF(nint(1.e6*gp_A))//"_mean"//COOP_STR_OF(nint(gp_mean))//"_width"//COOP_STR_OF(nint(gp_width))  , fwhm_arcmin = fwhm_arcmin)
#endif

contains

  !!define the fnl function of zeta; my convention of \zeta is:  \zeta = - \delta N (this determines the relative sign between \zeta and \Phi when I set up initial conditions in the Boltzmann code).

  function fnl_raw(x) result(y)
    COOP_SINGLE x, y
#ifdef BILLIARDS
    y =  deltaN%eval(abs(mean_chi+x*sigma_chi))  
#else
    y = gp_A* exp(-(x - gp_mean)**2/(2.d0*gp_width**2)) 
#endif
    return
  end function fnl_raw

  subroutine fnl(x, pixsize) 
    COOP_INT::n
    COOP_SINGLE x, y
    COOP_SINGLE pixsize, rms
    COOP_INT::i
    x = x/sigma_zeta
    if(pixsize .le. resolution)then
       x = fnl_raw(x)
       return
    endif
    rms = sqrt(3.*log(pixsize/resolution))
    n = max(2, ceiling(rms*2.5/0.1))
    y = fnl_raw(x)
    rms = rms*2.5/n
    do i = 1, n
       y = y + fnl_raw(x+i*rms)+fnl_raw(x-i*rms)
    enddo
    x = y/(2.*n+1.)
  end subroutine fnl

end program shells
