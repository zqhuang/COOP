program test
  use coop_wrapper_firstorder
  use coop_forecast_mod
  implicit none
#include "constants.h"
  type(coop_dataset_SN)::SN
  type(coop_MCMC_params)::mcmc

  call mcmc%init(prefix = "chains/wcdm", paramnames = "paramnames/wcdm.paramnames", ini = "myinis/wcdm.ini")
contains

  function loglike(mcmc)
    COOP_REAL loglike
    type(coop_MCMC_params)::mcmc
    loglike = SN%loglike(mcmc)
  end function loglike
  
  
  
end program test
