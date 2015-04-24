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
  COOP_STRING::inifile
  COOP_UNKNOWN_STRING, parameter::planckdata_path = "../data/cmb/"  ! "/home/zqhuang/includes/planck13/data" ! 
  COOP_INT i
  COOP_INT,parameter::total_steps = 6000
  COOP_INT,parameter::update_freq = 300
  logical do_update_propose 
  COOP_REAL::loglike
  call coop_MPI_init()

  if(iargc().ge.1)then
     inifile = trim(coop_InputArgs(1))
  else
     write(*,*) "missing ini file"
     write(*,*) "Syntax: "
     write(*,*) "./DOCLIK myinis/xx.ini"
     stop 
  endif

  mcmc%cosmology => cosmology
  call mcmc%init(trim(inifile))
  mcmc%do_flush = .true.  !!do not buffer
  
  if(mcmc%do_fastslow .and. coop_MPI_Rank().eq.0)then
     write(*,*) "doing fast-slow MCMC"
     write(*,*) "number of varying fast parameters:"//COOP_STR_OF(mcmc%n_fast)
     write(*,*) "number of varying slow parameters:"//COOP_STR_OF(mcmc%n_slow)
  endif

  
  call pl(1)%init(trim(planckdata_path)//"/CAMspec_v6.2TN_2013_02_26_dist.clik")
  call pl(2)%init(trim(planckdata_path)//"/commander_v4.1_lm49.clik")
  call pl(3)%init(trim(planckdata_path)//"/lowlike_v222.clik")  
  pool%CMB%cliklike => pl
  pool%HST%HSTlike => HSTlike


  call jla%read("../data/jla/jla.dataset")
  pool%SN_JLA%JLALike => jla
  do_update_propose = .true.
  
  !!do MCMC
  do i = 1, total_steps
     print*, "on Node ", coop_MPI_Rank(), ": step", i, " likelihood = ", mcmc%loglike
     if(do_update_propose)then
        if(mod(i, update_freq).eq.0)then
           call mcmc%update_propose()
        endif
        if(i .gt. total_steps/3)then
           do_update_propose = .false.
           call mcmc%chain%init()
           mcmc%do_memsave = .false.
        endif
     endif
     call mcmc%mcmc_step(pool)
  enddo
  call coop_MPI_finalize()
end program test
