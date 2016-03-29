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
  COOP_REAL::sigma_chi = 1.d-7
  COOP_REAL::mean_chi = 5.8d-6
  COOP_REAL::gp_width = 3.
  COOP_REAL::gp_A = 4.e-4
  COOP_REAL::gp_mean = 12.
  COOP_INT:: i, j

  call coop_file_load_function("deltaN-LUT-1.875", 1, 2, deltaN, .false.)
  call cosmology%Set_Planck_bestfit()
  call cosmology%compute_source(0)
  call coop_feedback("source done")
#ifdef BILLIARDS
  do i= 0 , 60, 5
     mean_chi = i*1.d-7
     write(*,*) "<chi> = ", mean_chi
     call coop_zeta3d_generate_cmb( cosmology, fnl, 300, 256, "zetaproj/testmap", "zetaproj/gp_meanchi"//COOP_STR_OF(nint(mean_chi*1.e7))//"_sigmachi"//COOP_STR_OF(nint(sigma_chi*1.e7)) )
  enddo
#else
!!$  do i=3, 3
!!$     print*, i
!!$     gp_width = i
  call coop_zeta3d_generate_cmb( cosmology, fnl, 300, 256, "zetaproj/testmap", "zetaproj/gs_amp"//COOP_STR_OF(nint(1.e6*gp_A))//"_mean"//COOP_STR_OF(nint(gp_mean))//"_width"//COOP_STR_OF(nint(gp_width)) )
!!$  enddo
#endif

contains

  !!define the fnl function of zeta; my convention of \zeta is:  \zeta = - \delta N (this determines the relative sign between \zeta and \Phi when I set up initial conditions in the Boltzmann code).

  subroutine fnl(x)
    real x
#ifdef BILLIARDS
    x = - deltaN%eval(abs(mean_chi+x/sigma_zeta*sigma_chi))  
#else
    x = - gp_A* exp(-(x/sigma_zeta - gp_mean)**2/(2.d0*gp_width**2)) 
#endif
    return
  end subroutine fnl

end program shells
