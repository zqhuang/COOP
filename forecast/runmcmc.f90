program RunMC
  use coop_wrapper_firstorder
  use coop_mcmc_mod
  implicit none
#include "constants.h"
  logical::use_CMB, use_SN, use_BAO, use_HST, use_WL, use_lensing, use_compressed_CMB, use_Age_Constraint
  type(coop_clik_object),target::pl(5)
  type(coop_HST_object),target::HSTlike
  type(coop_data_JLA),target:: jla
  type(coop_bao_object),target::bao(4)
  type(coop_wl_object), target::wl(1)
  type(coop_cosmology_firstorder),target::cosmology
  type(coop_dataset_CMB_simple),target::Compressed_CMB
  type(coop_dataset_Age_Constraint),target::Age
  type(coop_mcmc_params)::mcmc
  type(coop_data_pool)::pool
  type(coop_file)::fp
  type(coop_list_double)::z_out
  COOP_STRING::inifile, pname
  COOP_INT i, l, icmb, ik, iz
  COOP_REAL::loglike
  COOP_REAL::pvalue, norm, lnorm
  COOP_STRING::cls_root = ""
  COOP_STRING::cmb_dataset  = ""
  COOP_INT, parameter::nk = 256  
  COOP_REAL k(nk), matterPk(nk), khMpc(nk), WeylPk(nk)
  
  call coop_MPI_init()

  if(iargc().ge.1)then
     inifile = trim(coop_InputArgs(1))
  else
     write(*,*) "missing ini file"
     write(*,*) "Syntax: "
     write(*,*) "./MCMC myinis/xx.ini"
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
     call coop_dictionary_lookup(mcmc%settings, "use_CMB", use_CMB, .false.)
     call coop_dictionary_lookup(mcmc%settings, "use_WL", use_WL, .false.)     
     call coop_dictionary_lookup(mcmc%settings, "use_lensing", use_lensing, .false.)  
     call coop_dictionary_lookup(mcmc%settings, "use_SN", use_SN, .false.)  
     call coop_dictionary_lookup(mcmc%settings, "use_BAO", use_BAO, .false.)

     !!HST
     call coop_dictionary_lookup(mcmc%settings, "use_HST", use_HST, .false.)
     if(use_HST)then
        call coop_dictionary_lookup(mcmc%settings, "H0_center", HSTLike%H0, 70.6d0)
        call coop_dictionary_lookup(mcmc%settings, "H0_error", HSTLike%H0_err, 3.3d0)
     endif
     
     !!CMB derived background constraint (only used when use_CMB =  .false.)
     call coop_dictionary_lookup(mcmc%settings, "use_compressed_CMB", use_compressed_CMB, .false.)
     if(use_compressed_cmb .and. .not. use_CMB)then
        pool%CMB_Simple => compressed_CMB
        if(mcmc%feedback.gt.0)write(*,*) "Using CMB derived background constraint"        
     endif
     !!Age constraint
     call coop_dictionary_lookup(mcmc%settings, "use_age_constraint", use_Age_constraint, .false.)     
     if(use_Age_constraint)then
        pool%Age_Constraint => Age
        if(mcmc%feedback.gt.0)write(*,*) "Using age constraint"
     endif
     
     !!BAO
     if(use_BAO)then
        call bao(1)%init("%DATASETDIR%bao/sdss_6DF_bao.dataset")
        call bao(2)%init("%DATASETDIR%bao/sdss_MGS_bao.dataset")
        call bao(3)%init("%DATASETDIR%bao/sdss_DR11LOWZ_bao.dataset")
        call bao(4)%init("%DATASETDIR%bao/sdss_DR11CMASS_bao.dataset")
        if(mcmc%feedback.gt.0)write(*,*) "Using BAO"        
        pool%BAO%baolike => bao
        if(mcmc%init_level .lt. coop_init_level_set_background) mcmc%init_level = coop_init_level_set_background
     endif

     !!HST
     if(use_HST)then
        pool%HST%HSTlike => HSTlike
        if(mcmc%feedback.gt.0)write(*,*) "Using HST"
        if(mcmc%init_level .lt. coop_init_level_set_background) mcmc%init_level = coop_init_level_set_background
     endif

     !!supernova  
     if(use_SN)then
        if(mcmc%feedback.gt.0)write(*,*) "Using JLA"
        call jla%read("%DATASETDIR%jla/jla.dataset")
        pool%SN_JLA%JLALike => jla
        if(mcmc%init_level .lt. coop_init_level_set_background) mcmc%init_level = coop_init_level_set_background
     endif

     !!wl
     if(use_WL)then
        if(mcmc%feedback .gt. 0) write(*,*) "Using weak lensing data"
        call wl(1)%init("%DATASETDIR%weaklensing/CFHTLENS_6bin.dataset")
        pool%WL%WLLike => wl
        if(mcmc%init_level .lt. coop_init_level_set_pert) mcmc%init_level = coop_init_level_set_pert
     endif
     
     if(use_CMB)then
        icmb = 1
        do
           if(icmb .ge. size(pl)) stop "too many cmb datasets"
           call coop_dictionary_lookup(mcmc%settings, "cmb_dataset"//COOP_STR_OF(icmb), cmb_dataset)
           if(trim(cmb_dataset).eq."")then
              exit
           endif
           call pl(icmb)%init(trim(cmb_dataset))
           write(*, "(A)") "loading "//trim(cmb_dataset)
           icmb = icmb + 1              
        enddo
        if(use_lensing)then
           call coop_dictionary_lookup(mcmc%settings, "cmb_lensing_dataset", cmb_dataset)
           if(trim(cmb_dataset).ne."")then
              call pl(icmb)%init(trim(cmb_dataset))              
              write(*, "(A)") "loading "//trim(cmb_dataset)
           else
              stop "You need to set cmb_lensing_dataset for CMB lensing"
           endif
        else
           icmb = icmb - 1
        endif
        if(icmb .gt. 0)then
           if(mcmc%init_level .lt. coop_init_level_set_Cls) mcmc%init_level = coop_init_level_set_Cls           
           pool%CMB%cliklike => pl(1:icmb)
        endif
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
        if(mcmc%cosmology%h() .eq. 0.d0) then
           write(*,*) "-ln(likelihood) = \infty"
           stop "Model ruled out"
        else
           write(*,*) "h = ", mcmc%cosmology%h()
           write(*,"(A, F10.2, A)") "Age = ", mcmc%cosmology%AgeGyr(), " Gyr"           
           write(*,*) "omega_m = ", mcmc%cosmology%omega_m
           write(*,*) "sigma_8 = ", mcmc%cosmology%sigma_8
           write(*,*) "sigma_8 [rho_m/(3H^2)/0.3]^{0.3} = ", mcmc%cosmology%sigma_8 * (mcmc%cosmology%Omega_m*mcmc%cosmology%Mpsq0/0.3)**0.3           
           write(*,*) "omega_b h^2 M^2 = ", mcmc%cosmology%ombm2h2 
           write(*,*) "omega_c h^2 M^2 = ", mcmc%cosmology%omcm2h2
           write(*,*) "100theta = ", mcmc%cosmology%cosmomc_theta()*100.d0
           write(*,*) "z_recomb = ", mcmc%cosmology%zrecomb           
           write(*,*) "D_recomb = ", mcmc%cosmology%distlss
           
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
        if(iargc().ge.4)then
           call coop_get_Input(4, cls_root)
        else
           call coop_dictionary_lookup(mcmc%settings, "cls_root", cls_root)
        endif
        
        if(trim(cls_root) .ne. "" .and. use_CMB)then
           norm = mcmc%cosmology%TCmb()**2*1.d12           
           write(*,*) "saving Cl's to file: "//trim(cls_root)
           call fp%open(trim(cls_root)//"_scalCls.txt", "w")
           do l = 2, mcmc%lmax
              lnorm = l*(l+1.d0)/coop_2pi*norm
              write(fp%unit, "(I8, 20E16.7)") l, mcmc%cls_scalar(coop_index_ClTT, l)*lnorm, mcmc%cls_scalar(coop_index_ClEE, l)*lnorm, mcmc%cls_scalar(coop_index_ClTE, l)*lnorm, mcmc%cls_scalar(coop_index_ClLenLen, l)*norm*(l*(l+1.d0))**2, mcmc%cls_scalar(coop_index_ClTLen, l)*norm*(l*(l+1.d0))**1.5
           enddo
           call fp%close()
           if(mcmc%cosmology%has_tensor)then
              call fp%open(trim(cls_root)//"_tensCls.txt", "w")
              do l = 2,mcmc%lmax
                 lnorm = l*(l+1.d0)/coop_2pi*norm                 
                 write(fp%unit, "(I8, 20E16.7)") l, mcmc%cls_tensor(1:4, l)*lnorm
              enddo
              call fp%close()
           endif
           call fp%open(trim(cls_root)//"_lensedCls.txt", "w")
           do l = 2, mcmc%lmax
              lnorm = l*(l+1.d0)/coop_2pi             
              write(fp%unit, "(I8, 20E16.7)") l, mcmc%cls_lensed(1:4, l)*lnorm
           enddo
           call fp%close()
           call coop_set_uniform(nk, k, 0.4d0, 2.d3, logscale = .true.)
           khMpc = k * mcmc%cosmology%H0Mpc()/mcmc%cosmology%h()  !!k/H0 * (H0 * Mpc) / h = k in unit of h Mpc^{-1}
           call coop_dictionary_lookup(mcmc%settings, "output_redshifts", z_out)           
           do iz = 1, z_out%n
              !!compute k^3 |\delta_k|^2 /(2pi^2)  at redshift zero
              call mcmc%cosmology%get_Matter_power(z=z_out%element(iz), nk = nk, k = k, Pk = matterPk)
              matterPk = matterPk * (2.d0*coop_pi**2)/khMpc**3/mcmc%cosmology%h()**3  !!obtain |\delta_k|^2 in unit of (Mpc)^3
              call fp%open(trim(cls_root)//"_MatterPower"//COOP_STR_OF(iz)//".txt", "w")
              write(fp%unit, "(A, F10.2)") "#matter power at z = ", z_out%element(iz)
              do ik=1, nk
                 write(fp%unit, "(4E16.7)") khMpc(ik), matterPk(ik)
              enddo
           
              call fp%close()
           enddo
        endif
     else
        loglike = pool%loglike(mcmc)
        write(*,*) "-ln(likelihood) = ", loglike
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
