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
  type(coop_cosmology_firstorder)::cosmology
  type(coop_function)::deltaN
  call coop_file_load_function("deltaN-LUT-1.875", 1, 2, deltaN, .false.)
  call cosmology%Set_Planck_bestfit()
  call cosmology%compute_source(0)
  call coop_feedback("source done")
  call coop_3d_generate_cmb(cosmology, fnl, 300, 256, "zetaproj/testmap", "zetaproj/gp")

contains

  subroutine fnl(x)
    real x
    COOP_REAL,parameter::sigma_zeta = 4.6e-5
    COOP_REAL,parameter::sigma_chi = 5.e-7
    COOP_REAL,parameter::mean_chi = 5.8e-6
    x = deltaN%eval(abs(mean_chi+x/sigma_zeta*sigma_chi))
    return
  end subroutine fnl

end program shells
