program test
  use coop_wrapper_firstorder
  use coop_forecast_mod
  implicit none
#include "constants.h"
  type(coop_clik_object)::pl
  call coop_MPI_init()
  call pl%init("../data/cmb/hi_l/plik/plik_dx11dr2_HM_v18_TTTEEE.clik/")
  call pl%print_names()
  call pl%init("../data/cmb/low_l/bflike/lowl_SMW_70_dx11d_2014_10_03_v5c_Ap.clik/")
  call pl%print_names()
  call pl%init("../data/cmb/low_l/commander/commander_rc2_v1.1_l2_29_B.clik/")
  call pl%print_names()
  call pl%init("../data/cmb/lensing/smica_g30_ftl_full_pp.clik_lensing/")
  call pl%print_names()
  
  call coop_MPI_finalize()
end program test
