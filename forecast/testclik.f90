program test
  use coop_wrapper_firstorder
  use coop_forecast_mod
  implicit none
#include "constants.h"
  COOP_STRING::action= "TEST"
  logical::use_CMB, use_SN, use_BAO, use_HST, use_lensing, use_compressed_CMB
  type(coop_clik_object),target::pl(4)
  type(coop_HST_object),target::HSTlike
  type(coop_data_JLA),target:: jla
  type(coop_bao_object),target::bao(4)
  type(coop_cosmology_firstorder),target::cosmology
  type(coop_dataset_CMB_simple),target::Compressed_CMB
  type(coop_mcmc_params)::mcmc
  type(coop_data_pool)::pool
  COOP_STRING::inifile
  COOP_UNKNOWN_STRING, parameter::planckdata_path = "../data/cmb/"
  COOP_UNKNOWN_STRING, parameter::planckdata_path2 =   "/home/zqhuang/includes/planck13/data" 
  COOP_INT i
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
     if(mcmc%n_slow .eq.0)then
        mcmc%fast_per_round = 100
     endif
  endif

  
  call coop_dictionary_lookup(mcmc%settings, "action", action)
  if(trim(action) .eq. "") action = "TEST"  
  if(.not. mcmc%do_general_loglike)then
     call coop_dictionary_lookup(mcmc%settings, "use_CMB", use_CMB, .true.)
     call coop_dictionary_lookup(mcmc%settings, "use_lensing", use_lensing, .true.)  
     call coop_dictionary_lookup(mcmc%settings, "use_SN", use_SN, .true.)  
     call coop_dictionary_lookup(mcmc%settings, "use_BAO", use_BAO, .true.)
     call coop_dictionary_lookup(mcmc%settings, "use_HST", use_HST, .true.)
     call coop_dictionary_lookup(mcmc%settings, "use_compressed_CMB", use_compressed_CMB, .false.)



     if(use_compressed_cmb .and. .not. use_CMB)then
        pool%CMB_Simple => compressed_CMB
     endif

     !!BAO
     if(use_BAO)then
        call bao(1)%init("../data/bao/sdss_6DF_bao.dataset")
        call bao(2)%init("../data/bao/sdss_MGS_bao.dataset")
        call bao(3)%init("../data/bao/sdss_DR11LOWZ_bao.dataset")
        call bao(4)%init("../data/bao/sdss_DR11CMASS_bao.dataset")     
        pool%BAO%baolike => bao
     endif

     !!HST
     if(use_HST) pool%HST%HSTlike => HSTlike

     !!supernova  
     if(use_SN)then
        call jla%read("../data/jla/jla.dataset")
        pool%SN_JLA%JLALike => jla
     endif

     if(use_CMB)then
        if(coop_file_exists(trim(planckdata_path)//"/CAMspec_v6.2TN_2013_02_26_dist.clik"))then
           call pl(1)%init(trim(planckdata_path)//"/CAMspec_v6.2TN_2013_02_26_dist.clik")
           call pl(2)%init(trim(planckdata_path)//"/commander_v4.1_lm49.clik")
           call pl(3)%init(trim(planckdata_path)//"/lowlike_v222.clik")
           if(use_lensing)call pl(4)%init(trim(planckdata_path)//"/lensing_likelihood_v4_ref.clik_lensing")
        else !!try another path
           call pl(1)%init(trim(planckdata_path2)//"/CAMspec_v6.2TN_2013_02_26_dist.clik")
           call pl(2)%init(trim(planckdata_path2)//"/commander_v4.1_lm49.clik")
           call pl(3)%init(trim(planckdata_path2)//"/lowlike_v222.clik")
           if(use_lensing)call pl(4)%init(trim(planckdata_path2)//"/lensing_likelihood_v4_ref.clik_lensing")        

        endif
        pool%CMB%cliklike => pl
     endif
  endif
  select case(trim(action))
  case("TEST", "test")
     mcmc%params = mcmc%center
     mcmc%fullparams(mcmc%used) = mcmc%params
     if(associated(mcmc%cosmology))then
        call mcmc%get_lmax_from_data(pool)
        call mcmc%set_cosmology()
     endif        
     loglike = pool%loglike(mcmc)
     write(*,*) "-ln(likelihood) = ", loglike
  case("MCMC", "mcmc")
     if(mcmc%feedback .gt. 2) write(*,*) "Starting MCMC on Node #"//COOP_STR_OF(mcmc%proc_id)
     do i = 1, mcmc%total_steps
        call mcmc%mcmc_step(pool)
     enddo
  case default
     print*, trim(action)
     stop "unknown action"
  end select
  call coop_MPI_finalize()
end program test
