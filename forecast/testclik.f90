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
  COOP_UNKNOWN_STRING, parameter::planckdata_path = "/home/zqhuang/includes/planck13/data" ! "../data/cmb/"  !
  COOP_INT i
  COOP_INT,parameter::total_steps = 6000
  COOP_INT,parameter::update_freq = 300

  COOP_REAL::loglike
  call coop_MPI_init()
  
  mcmc%cosmology => cosmology
  call mcmc%init(prefix="chains/planck", paramnames= "paramnames/qcdm.paramnames", ini = "myinis/qcdm.ini")

  
  call pl(1)%init(trim(planckdata_path)//"/CAMspec_v6.2TN_2013_02_26_dist.clik")
  call pl(2)%init(trim(planckdata_path)//"/commander_v4.1_lm49.clik")
  call pl(3)%init(trim(planckdata_path)//"/lowlike_v222.clik")  
  pool%CMB%cliklike => pl
  pool%HST%HSTlike => HSTlike


  call jla%read("../data/jla/jla.dataset")
  pool%SN_JLA%JLALike => jla

!!$  call mcmc%set_cosmology()  
!!$  loglike = pool%loglike(mcmc)
!!$ print*, loglike
  !!do MCMC
  do i = 1, total_steps
     print*, "on Node ", coop_MPI_Rank(), ": step", i, " likelihood = ", mcmc%loglike
     if(i.lt. total_steps/4 .and. mod(i - update_freq/2, update_freq).eq.0)then
        call mcmc%update_propose()
        if(i.lt.update_freq)call mcmc%chain%init()
     endif
     call mcmc%mcmc_step(pool)
  enddo
  call coop_MPI_finalize()
end program test
