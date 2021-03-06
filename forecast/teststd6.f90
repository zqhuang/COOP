program test
  use coop_wrapper_firstorder
  use coop_forecast_mod
  implicit none
#include "constants.h"
  type(coop_clik_object),target::pl(1)
  type(coop_HST_object),target::HSTlike
  type(coop_data_JLA),target:: jla
  type(coop_cosmology_firstorder),target::cosmology

  type(coop_mcmc_params)::mcmc
  type(coop_data_pool)::pool
  COOP_UNKNOWN_STRING, parameter::planckdata_path = "../data/cmb/" !"/home/zqhuang/includes/planck13/data" !
  COOP_INT i

  COOP_REAL::loglike
  call coop_MPI_Init()
  mcmc%cosmology => cosmology
  call mcmc%init( "myinis/teststd6.ini")
  call mcmc%set_cosmology()
 call pl(1)%init(trim(planckdata_path)//"/CAMspec_v6.2TN_2013_02_26_dist.clik")
 ! call pl(1)%init(trim(planckdata_path)//"/commander_v4.1_lm49.clik")
!  call pl(1)%init(trim(planckdata_path)//"/lowlike_v222.clik")  
  pool%CMB%cliklike => pl
  loglike = pool%loglike(mcmc)
  write(*,*) "CMB loglike = ", loglike
!!$  nullify(pool%CMB%cliklike)
!!$  pool%HST%HSTlike => HSTlike
!!$  loglike = pool%loglike(mcmc)
!!$  write(*,*) "HST loglike = ", loglike
!!$  nullify(pool%HST%HSTlike)  
!!$
!!$  call jla%read("../data/jla/jla.dataset")
!!$  pool%SN_JLA%JLALike => jla
!!$  loglike = pool%loglike(mcmc)
!!$  write(*,*) "JLA loglike = ", loglike
  call coop_MPI_Finalize()
end program test
