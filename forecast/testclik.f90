program test
  use coop_wrapper_firstorder
  use coop_forecast_mod
  implicit none
#include "constants.h"
  type(coop_clik_object),target::pl(3)
  type(coop_HST_object),target::HSTlike
  type(coop_data_JLA),target:: jla
  type(coop_cosmology_firstorder),target::cosmology

  type(coop_mcmc_params)::mcmc
  type(coop_data_pool)::pool

  COOP_REAL::loglike
  mcmc%cosmology => cosmology
  call mcmc%init(prefix="chains/std6", paramnames= "paramnames/mcmc.paramnames", ini = "myinis/qcdm.ini")

  
  call pl(1)%init("../data/cmb/CAMspec_v6.2TN_2013_02_26_dist.clik")
  call pl(2)%init("../data/cmb/commander_v4.1_lm49.clik")
  call pl(3)%init("../data/cmb/lowlike_v222.clik")  
  pool%CMB%cliklike => pl
  pool%HST%HSTlike => HSTlike


  call jla%read("../data/jla/jla.dataset")
  pool%SN_JLA%JLALike => jla

  call mcmc%set_cosmology()  
  loglike = pool%loglike(mcmc)
  print*, loglike

end program test
