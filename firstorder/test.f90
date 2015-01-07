program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::fod
  type(coop_file)fp
  COOP_INT, parameter::nk=200
  COOP_INT ik
  COOP_REAL k(nk), Pk(nk), phi(nk), lnk(nk)
  call coop_set_uniform(nk, lnk, log(0.3d0), log(2.d3))
  k = exp(lnk)
  call fod%Set_Planck_bestfit()
  call fod%compute_source(0)
  print*, "sigma_8=",fod%sigma_8, fod%sigma_TopHat_R(0.d0, 8.d0*(1.d5/coop_SI_c)), fod%sigma_Gaussian_R(0.d0, 8.d0*(1.d5/coop_SI_c)), fod%sigma_Gaussian_R_quick(0.d0, 8.d0*(1.d5/coop_SI_c))
end program test
