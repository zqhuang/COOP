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
  call cosmology%Set_Planck_bestfit()
  call cosmology%compute_source(0)
  call coop_feedback("source done")
  call coop_3d_generate_cmb(cosmology, fnl, 10, 16, "testmap")

contains

  subroutine fnl(x)
    real x
    return
  end subroutine fnl

end program shells
