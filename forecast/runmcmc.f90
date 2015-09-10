program RunMC
  use coop_wrapper_firstorder
  use coop_forecast_mod
  implicit none
#include "constants.h"
  logical,parameter::use_Planck13 = .false.
  logical::use_CMB, use_SN, use_BAO, use_HST, use_lensing, use_compressed_CMB
  type(coop_clik_object),target::pl(4)
  type(coop_HST_object),target::HSTlike
  type(coop_data_JLA),target:: jla
  type(coop_bao_object),target::bao(4)
  type(coop_cosmology_firstorder),target::cosmology
  type(coop_dataset_CMB_simple),target::Compressed_CMB
  type(coop_mcmc_params)::mcmc
  type(coop_data_pool)::pool
  type(coop_file)::fp
  COOP_STRING::inifile, pname
  COOP_UNKNOWN_STRING, parameter::planckdata_path = "../data/cmb/"
  COOP_INT i
  COOP_REAL::loglike
  COOP_REAL::pvalue
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
  if(mcmc%do_fastslow .and. coop_MPI_Rank().eq.0 .and. trim(mcmc%action) .eq. "MCMC")then
     write(*,*) "doing fast-slow MCMC"
     write(*,*) "number of varying fast parameters:"//COOP_STR_OF(mcmc%n_fast)
     write(*,*) "number of varying slow parameters:"//COOP_STR_OF(mcmc%n_slow)
     if(mcmc%n_slow .eq.0)then
        mcmc%fast_per_round = 100
     endif
  endif

  
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
        if(mcmc%feedback.gt.0)write(*,*) "Using BAO"        
        pool%BAO%baolike => bao
     endif

     !!HST
     if(use_HST)then
        pool%HST%HSTlike => HSTlike
        if(mcmc%feedback.gt.0)write(*,*) "Using HST"
     endif

     !!supernova  
     if(use_SN)then
        if(mcmc%feedback.gt.0)write(*,*) "Using JLA"
        call jla%read("../data/jla/jla.dataset")
        pool%SN_JLA%JLALike => jla
     endif

     if(use_CMB)then
        if(use_Planck13)then
           if(mcmc%feedback.gt.0)write(*,*) "Using Planck 13 likelihood"
           if(coop_file_exists(trim(planckdata_path)//"/CAMspec_v6.2TN_2013_02_26_dist.clik/_mdb"))then
              call pl(1)%init(trim(planckdata_path)//"/CAMspec_v6.2TN_2013_02_26_dist.clik")
              call pl(2)%init(trim(planckdata_path)//"/commander_v4.1_lm49.clik")
              call pl(3)%init(trim(planckdata_path)//"/lowlike_v222.clik")
              if(use_lensing)call pl(4)%init(trim(planckdata_path)//"/lensing_likelihood_v4_ref.clik_lensing")
           else 
              write(*,*) "you need to make planck likelihood symbolic links to "//trim(planckdata_path)
              stop
           endif
        else  !!use Planck15
           if(mcmc%feedback.gt.0)write(*,*) "Using Planck 15 likelihood"           
           if(coop_file_exists(trim(planckdata_path)//"/hi_l/plik/plik_dx11dr2_HM_v18_TTTEEE.clik/_mdb"))then
              call pl(1)%init(trim(planckdata_path)//"/hi_l/plik/plik_dx11dr2_HM_v18_TTTEEE.clik")              
              call pl(2)%init(trim(planckdata_path)//"low_l/commander/commander_rc2_v1.1_l2_29_B.clik")
              call pl(3)%init(trim(planckdata_path)//"low_l/bflike/lowl_SMW_70_dx11d_2014_10_03_v5c_Ap.clik")
              if(use_lensing)call pl(4)%init(trim(planckdata_path)//"/lensing/smica_g30_ftl_full_pp.clik_lensing")
           else !!try another path
              write(*,*) "you need to make planck likelihood symbolic links to "//trim(planckdata_path)
              stop
           endif
        endif
        pool%CMB%cliklike => pl
     endif
  endif
  select case(trim(mcmc%action))
  case("TEST", "test")
     write(*,*) "Running TEST:"
     mcmc%params = mcmc%center
     mcmc%fullparams(mcmc%used) = mcmc%params
     if(iargc() .ge. 3)then
        pname = trim(coop_inputArgs(2))
        call coop_get_Input(3, pvalue)
        i=mcmc%paramnames%index(trim(pname))
        if(i.ne.0)then
           mcmc%fullparams(i) = pvalue
        else
           write(*,*) trim(pname)//" is not found"
           stop
        endif
     endif
     if(associated(mcmc%cosmology))then
        call mcmc%get_lmax_from_data(pool)
        write(*,*) "Setting up cosmology"
        call coop_prtsystime(.true.)        
        call mcmc%set_cosmology()
        call coop_prtsystime()             
     endif
     if(mcmc%cosmology%h() .eq. 0.d0) then
        write(*,*) "-ln(likelihood) = \infty"
     else
        write(*,*) "Computing likelihood"     
        call coop_prtsystime(.true.)
        loglike = pool%loglike(mcmc)
        call coop_prtsystime()     
        write(*,*) "-ln(likelihood) = ", loglike
        if(iargc().ge.3)then
           call fp%open(trim(mcmc%prefix)//".log", "a")        
           write(fp%unit, "(2G16.7)") pvalue, loglike        
        else
           call fp%open(trim(mcmc%prefix)//".log", "w")        
           write(fp%unit, "(G16.7)") loglike
        endif
        call fp%close()
     endif
  case("MCMC", "mcmc")
     if(mcmc%feedback .gt. 2) write(*,*) "Starting MCMC on Node #"//COOP_STR_OF(mcmc%proc_id)
     do i = 1, mcmc%total_steps
        call mcmc%mcmc_step(pool)
     enddo
  case("BEST", "best")
     call mcmc%findbest(pool, 0.01d0)
     call fp%open(trim(mcmc%prefix)//".best", "w")
     write(fp%unit, "(A16, G16.7)") "best like = ", mcmc%bestlike
     do i=1,mcmc%fulln
        write(fp%unit, "(A16, G16.7)")mcmc%paramnames%val(i), mcmc%fullparams(i)
     enddo
  case default
     print*, trim(mcmc%action)
     stop "unknown action"
  end select
  call coop_MPI_finalize()
end program RunMC
