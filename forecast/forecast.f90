module coop_forecast_mod
  use coop_wrapper_firstorder
  use coop_clik_mod
  use coop_HSTlike_mod
  use coop_SNlike_JLA_mod
  use coop_bao_mod
  use coop_wl_mod
  implicit none
#include "constants.h"

#define MCMC_OMEGA_M mcmc%fullparams(mcmc%index_omegam)
#define MCMC_OMEGA_K mcmc%fullparams(mcmc%index_omegak)
#define MCMC_W    mcmc%fullparams(mcmc%index_de_w)
#define MCMC_WA    mcmc%fullparams(mcmc%index_de_wa)
#define MCMC_OMEGA_LAMBDA  (1.d0 - MCMC_OMEGA_M - MCMC_OMEGA_K)  

  COOP_INT, parameter::coop_n_derived_with_cosmology = 4
  COOP_INT, parameter::coop_n_derived_without_cosmology = 1
  
  type coop_DataSet
     COOP_STRING::name
     type(coop_arguments)::args
     logical::off = .false.
     COOP_INT:: n_nuis = 0
     COOP_STRING, dimension(:),allocatable::nuis_names
   contains
     procedure::LogLike => coop_dataset_loglike
  end type coop_DataSet

  type, extends(coop_dataset):: coop_dataset_CMB_simple
     COOP_REAL::z_star = 1089.
     COOP_REAL::ombh2_center = 0.02207
     COOP_REAL::ombh2_sigma = 0.0004
     COOP_REAL::omch2_center = 0.1196
     COOP_REAL::omch2_sigma = 0.004
     COOP_REAL::theta_center = 1.04132
     COOP_REAL::theta_sigma = 0.0007
     COOP_REAL::corr_ombh2_omch2 = -0.558
     COOP_REAL::corr_ombh2_theta = 0.4675
     COOP_REAL::corr_omch2_theta = -0.4562
     COOP_REAL::invcov(3,3)
     logical::has_invcov = .false.
   contains
     procedure::loglike =>coop_dataset_CMB_simple_loglike
  end type coop_dataset_CMB_simple

  type, extends(coop_DataSet):: coop_dataset_SN_Simple
     COOP_INT::n = 0
     COOP_REAL,dimension(:),allocatable::z, mu, dmu, invdmusq
     COOP_REAL::suminvdmusq
     COOP_REAL::pec_vel = 400.d0/3.d5
   contains
     procedure::LogLike =>coop_dataset_SN_Simple_loglike
     procedure::import => coop_dataset_SN_Simple_import
     procedure::simulate => coop_dataset_SN_Simple_simulate
     procedure::export => coop_dataset_SN_Simple_export   
  end type coop_dataset_SN_Simple


  type, extends(coop_dataset)::coop_dataset_BAO
     type(coop_bao_object),dimension(:),pointer::baolike => null()
   contains
     procedure::loglike => coop_dataset_BAO_loglike
  end type coop_dataset_BAO
  
  type, extends(coop_dataset)::coop_dataset_CMB
     type(coop_clik_object),dimension(:),pointer::cliklike => null()
   contains
     procedure::LogLike => coop_dataset_CMB_LogLike
  end type coop_dataset_CMB

  type, extends(coop_dataset)::coop_dataset_HST
     type(coop_HST_object), pointer::HSTlike => null()
   contains
     procedure::LogLike => coop_dataset_HST_logLike
  end type coop_dataset_HST


  type, extends(coop_dataset)::coop_dataset_SN_JLA
     type(coop_data_JLA), pointer::JLAlike => null()
   contains
     procedure::loglike => coop_dataset_SN_JLA_loglike
  end type coop_dataset_SN_JLA

  type, extends(coop_dataset)::coop_dataset_WL
     type(coop_WL_object), dimension(:), pointer::WLLike => null()
   contains
     procedure::LogLike => coop_dataset_WL_LogLike
  end type coop_dataset_WL


  type coop_Data_Pool
!!these are simulations     
     type(coop_dataset_SN_Simple),dimension(:), pointer::SN_Simple => null()
     type(coop_dataset_CMB_Simple), pointer::CMB_Simple => null()
     type(coop_dataset_BAO)::BAO
     type(coop_dataset_CMB)::CMB
     type(coop_dataset_HST)::HST
     type(coop_dataset_SN_JLA)::SN_JLA
     type(coop_dataset_WL)::WL
   contains
     procedure::LogLike => coop_data_pool_LogLike
  end type coop_Data_Pool

  
  type coop_MCMC_params
     COOP_STRING::action
     COOP_INT::feedback = 0
     COOP_INT::proc_id = 0
     COOP_INT::total_steps = 50000
     type(coop_cosmology_firstorder),pointer::cosmology => null()
     type(coop_cosmology_firstorder),pointer::cosmology_saved => null()        
     type(coop_file)::chainfile
     type(coop_file)::ndffile     
     type(coop_dictionary)::settings     
     COOP_STRING::prefix, chainname, ndfname
     COOP_STRING::form
     logical::do_flush = .false.
     logical::do_write_chain = .true.
     logical::do_ndf = .true.  !!for likelihood interpolation
     logical::do_drift = .false.
     logical::do_overwrite = .false.
     logical::do_fastslow = .false.
     logical::slow_changed = .true.
     logical::do_memsave = .true.
     logical::do_general_loglike = .false.
     COOP_REAL::converge_R = coop_logZero
     COOP_REAL::approx_frac = 0.d0
     COOP_REAL::drift_frac = 0.d0
     COOP_REAL::drift_step = 0.01d0
     logical::is_drift = .false.
     COOP_INT::n_fast = 0
     COOP_INT::index_fast_start = 0
     COOP_INT::n_slow = 0
     COOP_INT::index_propose_fast = 1
     COOP_INT::index_propose_slow = 1     
     COOP_INT::index_propose = 1
     COOP_INT::fast_per_round = 5
     COOP_INT::time = 0
     COOP_INT:: n = 0
     COOP_INT:: fulln = 0
     COOP_INT:: n_derived = 0
     COOP_REAL:: proposal_length = 2.4d0
     COOP_REAL::bestlike = coop_LogZero
     COOP_REAL::loglike = coop_LogZero
     logical::loglike_is_exact = .false.
     COOP_REAL::loglike_proposed = coop_LogZero
     logical::loglike_proposed_is_exact = .false.     
     COOP_REAL::mult = 0.d0
     COOP_REAL::sum_mult = 0.d0
     COOP_REAL::temperature = 1.d0
     COOP_INT::lmax = 0
     COOP_INT::update_seconds = 0 !!this has higher priority unless not used (set to be zero)
     COOP_INT::update_steps = 0  !!lower priority; not used if set to be zero
     COOP_INT::num_exact_calc = 0
     COOP_INT::num_approx_calc = 0     
     COOP_REAL, dimension(:,:), allocatable::Cls_scalar, Cls_tensor, Cls_lensed     
     COOP_INT,dimension(:),allocatable::used, map2used
     COOP_REAL,dimension(:),allocatable::fullparams
     COOP_REAL,dimension(:),allocatable::params
     COOP_REAL,dimension(:),allocatable::bestparams     
     COOP_SHORT_STRING, dimension(:), allocatable::name
     COOP_SHORT_STRING, dimension(:), allocatable::tex
     COOP_REAL,dimension(:),allocatable::lower
     COOP_REAL,dimension(:),allocatable::upper
     COOP_REAL,dimension(:),allocatable::center
     COOP_REAL,dimension(:),allocatable::width
     logical,dimension(:),allocatable::has_prior
     logical::any_prior = .false.
     COOP_INT::n_lc_priors = 0
     COOP_REAL,dimension(:),allocatable::prior_sigma, lc_prior_sigma     
     COOP_REAL,dimension(:),allocatable::prior_center, lc_prior_center
     COOP_INT,dimension(:,:), allocatable::lc_prior_index
     COOP_REAL,dimension(:,:), allocatable::lc_prior_weight     
     COOP_REAL,dimension(:,:),allocatable::vecs
     COOP_REAL,dimension(:,:),allocatable::vecs_fast
     COOP_REAL,dimension(:,:),allocatable::vecs_slow     
     COOP_REAL,dimension(:),allocatable::params_saved
     COOP_SINGLE, dimension(:),allocatable::knot
     COOP_REAL,dimension(:,:),allocatable::mapping     
     COOP_REAL,dimension(:,:),allocatable::mapping_fast
     COOP_REAL,dimension(:),allocatable::derived_params
     type(coop_random_cycl)::cycl
     type(coop_covmat)::covmat
     type(coop_list_realarr)::chain
     type(coop_dictionary)::paramnames
     type(coop_nd_prob)::like_approx
     COOP_INT::accept = 0
     COOP_INT::reject = 0
     COOP_INT::step = 0
     !!index for parameters
     COOP_INT::index_ombm2h2 = 0
     COOP_INT::index_omcm2h2 = 0          
     COOP_INT::index_theta = 0
     COOP_INT::index_tau = 0
     COOP_INT::index_mnu = 0
     COOP_INT::index_logA = 0
     COOP_INT::index_logAm2tau = 0     
     COOP_INT::index_ns = 0
     COOP_INT::index_nrun = 0     
     COOP_INT::index_r = 0
     COOP_INT::index_nt = 0
     COOP_INT::index_omegam = 0
     COOP_INT::index_omegab = 0     
     COOP_INT::index_omegak = 0
     COOP_INT::index_h = 0
     COOP_INT::index_de_w = 0
     COOP_INT::index_de_wa = 0
     COOP_INT::index_de_epss = 0
     COOP_INT::index_de_epsinf = 0
     COOP_INT::index_de_zetas = 0
     COOP_INT::index_de_betas = 0                         
#if DO_COUPLED_DE     
     COOP_INT::index_de_Q = 0
#elif DO_EFT_DE     
     COOP_INT::index_de_alpha_M0 = 0
     COOP_INT::index_de_alpha_H0 = 0
     COOP_INT::index_de_alpha_T0 = 0
     COOP_INT::index_de_alpha_B0 = 0
     COOP_INT::index_de_alpha_K0 = 0
     logical::w_is_background = .false.
#endif     
   contains
     procedure::init => coop_MCMC_params_init
     procedure::MCMC_step => coop_MCMC_params_MCMC_step
     procedure::update_propose => coop_MCMC_params_update_Propose
     procedure::findbest => coop_MCMC_params_findbest
     procedure::derived => coop_MCMC_params_derived
     procedure::priorLike => coop_MCMC_params_priorLike
     procedure::Set_Cosmology => coop_MCMC_params_Set_Cosmology
     procedure::index_of => coop_MCMC_params_index_of
     procedure::get_lmax_from_data => coop_MCMC_params_get_lmax_from_data
     procedure::proposal_r => coop_MCMC_params_proposal_r
     procedure::propose_vec => coop_MCMC_params_propose_vec
     procedure::propose_fast_vec => coop_MCMC_params_propose_fast_vec
     procedure::propose_slow_vec => coop_MCMC_params_propose_slow_vec     
     procedure::general_loglike => coop_MCMC_params_general_loglike
  end type coop_MCMC_params



contains

  function coop_MCMC_params_index_of(this, name) result(ind)
    class(coop_MCMC_params)::this    
    COOP_UNKNOWN_STRING::name
    COOP_INT::ind
    ind = this%paramnames%index(trim(adjustl(name)))
  end function coop_MCMC_params_index_of

  function coop_MCMC_params_priorLike(this) result(Prior)
    class(coop_MCMC_params)::this
    COOP_REAL::Prior
    COOP_INT::ip
    Prior = coop_LogZero
    if(any(this%params .gt. this%upper) .or. any(this%params.lt.this%lower))return
    if(this%any_prior)then
       Prior = sum(((this%params - this%prior_center)/this%prior_sigma)**2, mask = this%has_prior)/2.d0
    else
       prior = 0.d0
    endif
    if(this%n_lc_priors .gt. 0)then
       do ip = 1, this%n_lc_priors
          Prior = Prior + ((this%fullparams(this%lc_prior_index(1, ip)) * this%lc_prior_weight(1, ip) + this%fullparams(this%lc_prior_index(2, ip)) * this%lc_prior_weight(2, ip) - this%lc_prior_center(ip))/this%lc_prior_sigma(ip))**2/2.d0
       enddo
    endif
  end function coop_MCMC_params_priorLike

  !!set up %cosmology from %fullparams
  subroutine coop_MCMC_params_Set_Cosmology(this)
    class(coop_MCMC_params)::this
    COOP_REAL,parameter::omega_m_min = 0.15
    COOP_REAL,parameter::omega_m_max = 0.45   
    COOP_REAL, parameter::h_b_i = 0.45d0
    COOP_REAL, parameter::h_t_i =  0.95d0
    COOP_REAL, parameter::alpha_dep_Omega_m = 0.3d0
    COOP_REAL, parameter::alpha_dep_Omega_r = 8.d-5
    COOP_INT::iloop
    COOP_REAL::h_t, h_b, h_m, theta_t, theta_b, theta_m, theta_want
    logical success
    if(.not. associated(this%cosmology)) stop "MCMC_params_set_cosmology: cosmology not allocated"
#if DO_EFT_DE
    if(this%index_de_alpha_M0 .ne. 0)then    
       call this%cosmology%set_alphaM(coop_function_constructor( coop_de_alpha_invh2, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., &
            args= coop_arguments_constructor( r = (/ this%fullparams(this%index_de_alpha_M0), alpha_dep_Omega_m, alpha_dep_Omega_r /) ), name = "alpha_M" ))
    else
       call this%cosmology%set_alphaM( coop_function_polynomial( (/ 0.d0 /) ) )       
    endif
#endif

    if(this%index_theta .ne. 0)then
       theta_want = this%fullparams(this%index_theta)
       h_t = min(h_t_i, sqrt((this%fullparams(this%index_ombm2h2)+this%fullparams(this%index_omcm2h2))/omega_m_min/coop_Mpsq0))
       h_b = max(h_b_i, sqrt((this%fullparams(this%index_ombm2h2)+this%fullparams(this%index_omcm2h2))/omega_m_max/coop_Mpsq0))
       call calc_theta(h_t, theta_t)
       if(this%cosmology%h().eq.0.d0)return
       call calc_theta(h_b, theta_b)
       if(this%cosmology%h().eq.0.d0)return       
       if(theta_t .lt. theta_want)then
          call this%cosmology%set_h(0.d0)
          return
       endif
       if(theta_b .gt. theta_want)then
          call this%cosmology%set_h(0.d0)
          return
       endif
       iloop = 0
       do while(h_t - h_b .gt. 2.d-4)
          h_m = (h_t + h_b)/2.d0
          call calc_theta(h_m, theta_m)
          if(this%cosmology%h().eq.0.d0)return          
          if(theta_m .gt. theta_want)then
             h_t = h_m
             theta_t = theta_m
          else
             h_b = h_m
             theta_b = theta_m
          endif
          iloop = iloop+1
          if(iloop.gt. 20)then
             call this%cosmology%set_h(0.d0)
             if(this%feedback .ge. 1)write(*,*) "SetH on Node "//COOP_STR_OF(this%proc_id)//" failed"
             return
          endif
       enddo
       if((theta_t-theta_b).gt. 1.d-6)then
          h_m = ((theta_t - theta_want)*h_b + (theta_want - theta_b)*h_t)/(theta_t-theta_b)
          call setForH(h_m)
          if(this%cosmology%h().eq.0.d0)return
       endif
    elseif(this%index_h .ne. 0)then
       call setForH(this%fullparams(this%index_h))
          if(this%cosmology%h().eq.0.d0)return       
    else
       stop "you need to use either theta or h for MCMC runs"
    endif
#if DO_EFT_DE
    this%cosmology%f_alpha_M = coop_EFT_DE_alphaM    
    if(this%index_de_alpha_K0 .ne. 0)then
       this%cosmology%f_alpha_K = coop_function_constructor( coop_de_alpha_invh2, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args= coop_arguments_constructor ( r = (/ this%fullparams(this%index_de_alpha_K0), alpha_dep_Omega_m, this%cosmology%Omega_radiation() + this%cosmology%Omega_massless_neutrinos() /) ) , name = "EFT DE alpha_K")
    endif
    if(this%index_de_alpha_B0 .ne. 0)then
       this%cosmology%f_alpha_B = coop_function_constructor( coop_de_alpha_invh2, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args= coop_arguments_constructor ( r = (/ this%fullparams(this%index_de_alpha_B0), alpha_dep_Omega_m, this%cosmology%Omega_radiation() + this%cosmology%Omega_massless_neutrinos() /) ) , name = "EFT DE alpha_B")
    endif
    if(this%index_de_alpha_H0 .ne. 0)then
       this%cosmology%f_alpha_H = coop_function_constructor( coop_de_alpha_invh2, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args= coop_arguments_constructor ( r = (/ this%fullparams(this%index_de_alpha_H0), alpha_dep_Omega_m, this%cosmology%Omega_radiation() + this%cosmology%Omega_massless_neutrinos() /) ) , name = "EFT DE alpha_H")
    endif
    if(this%index_de_alpha_T0 .ne. 0)then
       this%cosmology%f_alpha_T = coop_function_constructor( coop_de_alpha_invh2, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args= coop_arguments_constructor ( r = (/ this%fullparams(this%index_de_alpha_T0), alpha_dep_Omega_m, this%cosmology%Omega_radiation() + this%cosmology%Omega_massless_neutrinos() /) ) , name = "EFT DE alpha_T")
    endif
#endif    

    call this%cosmology%setup_background()
    if(this%index_tau .ne. 0)then !!
       this%cosmology%optre = this%fullparams(this%index_tau)
       call this%cosmology%set_xe()
       if(this%index_logAm2tau .ne. 0)then
          this%cosmology%As = exp(this%fullparams(this%index_logAm2tau)+2.d0*this%cosmology%optre)*1.d-10      
       else
          if(this%index_logA .ne. 0)then
             this%cosmology%As = exp(this%fullparams(this%index_logA))*1.d-10
          else
             this%cosmology%As = 2.2d-9
          endif
       endif
       if(this%index_ns .ne. 0)then
          this%cosmology%ns = this%fullparams(this%index_ns)
       else
          this%cosmoLogy%ns = 0.967d0
       endif
       if(this%index_nrun .ne.0)then
          this%cosmology%nrun = this%fullparams(this%index_nrun)
       else
          this%cosmology%nrun = 0.d0
       endif
       if(this%index_r .ne. 0)then
          this%cosmology%r = this%fullparams(this%index_r)
          this%cosmology%has_tensor = (this%cosmology%r .gt. 0.d0)
       else
          this%cosmology%r = 0.d0
          this%cosmology%has_tensor = .false.
       endif
       if(this%index_nt .ne. 0)then
          this%cosmology%nt = this%fullparams(this%index_nt)
       else
          this%cosmology%nt = - this%cosmology%r/8.d0 !!inflationary consistency
       endif
       call this%cosmology%set_standard_power(this%cosmology%As, this%cosmology%ns, this%cosmology%nrun, this%cosmology%r, this%cosmology%nt)
       if(this%lmax .gt. 1)then
          call this%cosmology%compute_source(0, success = success)
          if(.not. success)then
             call this%cosmology%set_h(0.d0)
             return
          endif
          if(this%cosmology%has_tensor)then
             call this%cosmology%compute_source(2, success = success)
             if(.not. success)then
                call this%cosmology%set_h(0.d0)
                return
             endif
          endif
          if(this%feedback .gt. 2) then
             write(*,*) "sources done"
             call coop_prtsystime()
          endif
          
          call this%cosmology%source(0)%get_all_cls(2, this%lmax, this%Cls_scalar)
          call coop_get_lensing_Cls(2, this%lmax, this%Cls_Scalar, this%Cls_lensed)
          this%Cls_lensed = this%Cls_lensed + this%Cls_scalar
          if(this%cosmology%has_tensor)then
             call this%cosmology%source(2)%get_all_cls( 2, this%lmax, this%Cls_tensor)
             this%Cls_lensed = this%Cls_lensed + this%Cls_tensor
          endif
          this%Cls_lensed = this%Cls_lensed*((this%cosmology%Tcmb())**2*1.d12)
          if(this%feedback .gt. 2) then
             write(*,*) "lensing Cls done"
             call coop_prtsystime()
          endif
       endif       
    endif
    
  contains

    subroutine calc_theta(h, theta)
      COOP_REAL::h, theta
      call setForH(h)
      if(this%cosmology%h().eq.0.d0)return
      theta = this%cosmology%cosmomc_theta()*100.d0
    end subroutine calc_theta

    subroutine setforH(h)
      COOP_REAL::h
      COOP_REAL::Q, w, wa, epsilon_s, epsilon_inf, zeta_s, beta_s
      COOP_INT::err
      type(coop_function)::fQ, fwp1
      call this%cosmology%free()
      call this%cosmology%init(name = "Cosmology", id = 0, h = h)
      this%cosmology%ombh2 =this%fullparams(this%index_ombm2h2)/coop_Mpsq0
      this%cosmology%omch2 =this%fullparams(this%index_omcm2h2)/coop_Mpsq0
      !!baryon
      call this%cosmology%add_species(coop_baryon(this%cosmology%ombh2/h**2))
      this%cosmology%index_baryon = this%cosmology%num_species
      !!radiation
      call this%cosmology%add_species(coop_radiation(this%cosmology%Omega_radiation()))
      this%cosmology%index_radiation = this%cosmology%num_species      
      !!neutrinos
      if(this%index_mnu .ne. 0)then
#if DO_EFT_DE
         stop "For EFT DE model massive neutrinos are not implemented yet."
#endif         
         if(this%fullparams(this%index_mnu).gt. 0.01d0)then  !!do massive neutrinos
            call this%cosmology%add_species(coop_neutrinos_massive( &
                 this%cosmology%Omega_nu_per_species_from_mnu_eV( this%fullparams(this%index_mnu) ) ,&
                 this%cosmology%Omega_massless_neutrinos_per_species()))
            this%cosmology%index_massivenu = this%cosmology%num_species                  
            call this%cosmology%add_species( coop_neutrinos_massless(this%cosmology%Omega_massless_neutrinos_per_species()*(this%cosmology%NNu()-1)))
            this%cosmology%index_nu = this%cosmology%num_species                  
         else
            call this%cosmology%add_species( coop_neutrinos_massless(this%cosmology%Omega_massless_neutrinos()))
            this%cosmology%index_nu = this%cosmology%num_species            
         endif
      else
         call this%cosmology%add_species( coop_neutrinos_massless(this%cosmology%Omega_massless_neutrinos()))
         this%cosmology%index_nu = this%cosmology%num_species
      endif

      if(this%index_de_epss .ne. 0 .or. this%index_de_epsinf .ne. 0 .or. this%index_de_zetas .ne. 0 .or. this%index_de_betas .ne. 0)then
         if(this%index_de_epss .ne. 0)then
            epsilon_s = this%fullparams(this%index_de_epss)
         else
            epsilon_s = 0.d0
         endif
         if(this%index_de_epsinf .ne. 0)then
            epsilon_inf = this%fullparams(this%index_de_epsinf)
         else
            epsilon_inf = 0.01d0
         endif
         if(this%index_de_zetas .ne. 0)then
            zeta_s = this%fullparams(this%index_de_zetas)
         else
            zeta_s = 0.d0
         endif
         if(this%index_de_betas .ne. 0)then
            beta_s = this%fullparams(this%index_de_betas)
         else
            beta_s = 6.d0
         endif
         fwp1 = coop_function_constructor(coop_de_wp1_coupled_quintessence, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args = coop_arguments_constructor( r = (/ this%cosmology%Omega_k()  - this%cosmology%omch2/h**2, epsilon_s, epsilon_inf, zeta_s , beta_s /) ), name = "DE 1+w")
      else
         if(this%index_de_w .ne. 0)then
            w = this%fullparams(this%index_de_w)
            if(this%index_de_wa .ne. 0)then
               wa = this%fullparams(this%index_de_wa)
            else
               wa = 0.d0
            endif
         else
            w = -1.d0
            wa = 0.d0
         endif
         call fwp1%init_polynomial( (/ 1.d0 + w + wa, - wa /) )
      endif
      
#if DO_EFT_DE
      call this%cosmology%add_species(coop_cdm(this%cosmology%omch2/h**2))
      this%cosmology%index_cdm = this%cosmology%num_species      
      if(this%index_de_alpha_M0 .ne. 0)then
         if(this%w_is_background)then
            call coop_background_add_EFT_DE_with_effective_w(this%cosmology, effective_wp1 = fwp1 , err = err)            
         else
            call coop_background_add_EFT_DE(this%cosmology, wp1 = fwp1 , err = err)
         endif

         
         if(err .ne. 0) then
            call this%cosmology%set_h(0.d0)
            call fwp1%free()
            return
         endif
      else
         call this%cosmology%add_species(coop_de_general(this%cosmology%Omega_k(), fwp1))
      endif
         
      
      this%cosmology%de_genre = COOP_PERT_EFT
#elif DO_COUPLED_DE
      if(this%index_de_Q .ne. 0)then
         Q = this%fullparams(this%index_de_Q)
      else
         Q = 0.d0
      endif
      call fQ%init_polynomial( (/ Q /) )
      call coop_background_add_coupled_DE(this%cosmology, Omega_c = this%cosmology%omch2/h**2, fwp1 = fwp1, fQ = fQ, err = err)
      this%cosmology%index_cdm = this%cosmology%num_species-1                  
      if(err .ne. 0)then
         call this%cosmology%set_h(0.d0)
      endif
      this%cosmology%de_genre = COOP_PERT_SCALAR_FIELD      
      call fQ%free()
#else      
      call this%cosmology%add_species(coop_cdm(this%cosmology%omch2/h**2))
      this%cosmology%index_cdm = this%cosmology%num_species                  
      if(this%index_de_w .ne. 0)then
         w = this%fullparams(this%index_de_w)
         if(this%index_de_wa .ne. 0)then
            wa = this%fullparams(this%index_de_wa)
            call this%cosmology%add_species(coop_de_w0wa(this%cosmology%Omega_k(), w, wa))               
         else
            call this%cosmology%add_species(coop_de_w0(this%cosmology%Omega_k(), w))                              
         endif
      else
         call this%cosmology%add_species(coop_de_lambda(this%cosmology%Omega_k()))                           
      endif
#endif
      this%cosmology%index_de = this%cosmology%num_species                  
100   call fwp1%free()
    end subroutine setforH

  end subroutine coop_MCMC_params_Set_Cosmology

  function  coop_MCMC_params_derived(this) result(derived)
    class(coop_MCMC_params)::this
    COOP_REAL:: derived(this%n_derived), phi, dlnVdphi, aeq, Q
    if(this%n_derived .le. 0) return
    if(associated(this%cosmology))then
       derived(1) = this%cosmology%h()*100.d0
       derived(2) = this%cosmology%Omega_m
       derived(3) = 1.d0 - this%cosmology%Omega_m
       derived(4) = this%cosmology%sigma_8
    else
       if(this%index_omegam .ne. 0 .and. this%index_omegak .ne. 0)then
          derived(1) = 1.d0 - this%fullparams(this%index_omegam)- this%fullparams(this%index_omegak)
       else
          derived(1) = 0.d0
       endif
    endif
  end function coop_MCMC_params_derived

  subroutine coop_MCMC_params_update_Propose(this)
    class(coop_MCMC_params)::this
    COOP_INT::i, istart, i1, i2
    type(coop_file)::fp
    COOP_REAL:: mult, diff(this%n)
    if(this%step .eq. this%total_steps)then
       if(this%update_seconds .ne. 0) call coop_MPI_Abort("Step = Total Steps; MPI terminating.")
    endif
    if(this%update_seconds .gt. 0)then
       this%time = nint(coop_systime_sec())       
       if(this%time .lt. this%update_seconds)return
    else
       if(mod(this%step, this%update_steps).ne. 0)return
    endif
    if(this%feedback .ge. 1)write(*,*) "updating propose matrix on Node "//COOP_STR_OF(this%proc_id)//", step "//COOP_STR_OF(this%step)//", time = "//COOP_STR_OF(this%time)
    if(.not. this%do_memsave)stop "cannot update propose matrix when do_memsave is off"
    call this%covmat%alloc(this%n)  !!set sigma = 1 and the rest 0
    if(this%chain%n .ge. 10)then
       istart =  ceiling(this%chain%n*0.2)  !!discard 20% samples
       do i = istart, this%chain%n  
          call this%chain%get_element(i, this%knot)
          this%covmat%mult  = this%covmat%mult + this%knot(1)
          this%covmat%mean = this%covmat%mean + this%knot(1)*this%knot(3:this%n+2)
       enddo
       if(this%covmat%mult .gt. 1.d-2)then
          this%covmat%mean = this%covmat%mean / this%covmat%mult
          do i = istart, this%chain%n  
             call this%chain%get_element(i, this%knot)
             diff  = this%knot(3:this%n+2)  - this%covmat%mean
             do i1 = 1, this%n
                do i2 = 1, i1
                   this%covmat%c(i1, i2) = this%covmat%c(i1, i2) + diff(i1)*diff(i2)*this%knot(1)
                enddo
             enddo
          enddo
          do i1 = 1, this%n
             do i2 = 1, i1 -1         
                this%covmat%c(i2, i1) = this%covmat%c(i1, i2)
             enddo
          enddo
          this%covmat%c = this%covmat%c/this%covmat%mult
       endif
    endif
    if(this%do_ndf) call this%ndffile%close()
    call this%covmat%MPI_Sync(converge_R = this%converge_R) !!terminate the program       
    if(this%do_ndf)then
       call this%like_approx%load_chains(this%prefix, this%n)
       if(this%feedback .ge. 1)write(*,*) COOP_STR_OF(this%proc_id)//": likelihood fitting function loaded with "//COOP_STR_OF(this%like_approx%n)//" data points"
       call this%ndffile%open(this%ndfname, "ua")       
    endif
    this%time = nint(coop_systime_sec(0.d0))  !!reset time
    if(this%feedback.ge.1 .and. this%proc_id .eq. 0)then
       write(*, "(A, G15.4)") "convergence R - 1 = ", this%converge_R       
    endif
    if(this%proc_id.eq.0)then
       call fp%open(trim(this%prefix)//".converge_stat", "w")
       write(fp%unit, "(G15.4)") this%converge_R
       call fp%close()
    endif
       
       
    if(this%covmat%mult .gt. this%n*10.d0 .and. this%converge_R .lt. 1000.d0 .and. this%converge_R .gt. 0.03d0 .and. .not. coop_isnan(this%covmat%L) .and. all(this%covmat%sigma.gt.0.d0))then !!update mapping matrix
       if(this%proc_id.eq.0)call this%covmat%export(trim(this%prefix)//".runcov")       
       do i=1, this%n
          this%mapping(i, :) = this%covmat%L(i, :)*this%covmat%sigma(i)
       enddo
       if(this%do_fastslow)then
          this%mapping_fast = this%mapping(this%index_fast_start:this%n, this%index_fast_start:this%n)
       endif
       if(this%feedback .ge. 1)write(*,*) "propose matrix is updated on Node "//COOP_STR_OF(this%proc_id)
    endif
  end subroutine coop_MCMC_params_update_Propose


  

  subroutine coop_MCMC_params_MCMC_step(this, Pool)
    class(coop_MCMC_params)::this
    type(coop_data_pool)::pool
    COOP_REAL vec(this%n)
    COOP_INT i
    COOP_LONG_STRING::line
    logical::cosmology_changed
    type(coop_cosmology_firstorder),target::cosmo
    if(this%n .le. 0) call coop_MPI_Abort("MCMC: no varying parameters")
    select type(this)
    type is(coop_MCMC_params)
       if(this%feedback .gt. 4)then !!check reject rate
          if(this%reject .gt. this%accept*30 .and. this%accept .gt. 0)then
             write(*,*) "Reject rate enormalously high"
             write(*,*) "MPI Rank = ",this%proc_id, size(this%covmat%mean), size(this%covmat%sigma), size(this%covmat%c), size(this%covmat%L)             
             write(*,*) "mean = "             
             write(*,"("//COOP_STR_OF(this%n)//"E14.5)") this%covmat%mean
             write(*,*) "sigma = "                          
             write(*,"("//COOP_STR_OF(this%n)//"E14.5)") this%covmat%sigma
             write(*,*) "correlation matrix = "
             do i=1, this%covmat%n
                write(*,"("//COOP_STR_OF(this%n)//"E14.5)") this%covmat%c(i, :)
             enddo             
             write(*,*) "mapping matrix = "
             do i=1, this%covmat%n
                write(*,"("//COOP_STR_OF(this%n)//"E14.5)") this%mapping(i, :)
             enddo
          endif
       endif
       if(this%step .eq. 0)then
          call coop_random_init()
          call this%chain%init()
          this%sum_mult = 0.d0
          call this%get_lmax_from_data(pool)
          this%chainname = trim(this%prefix)//"_"//COOP_STR_OF(this%proc_id+1)//".txt"
          this%ndfname = trim(this%prefix)//"_"//COOP_STR_OF(this%proc_id+1)//".ndf"          
          if(this%do_overwrite)then
             if(this%do_ndf) call this%ndffile%open(this%ndfname, "u")
             
          else
             if(this%do_ndf)then
                call this%like_approx%load_chains(this%prefix, this%n)
                if(this%feedback .ge.2 .and. this%like_approx%n .gt. 0)then
                   write(*,*) COOP_STR_OF(this%proc_id)//": likelihood fitting function loaded with "//COOP_STR_OF(this%like_approx%n)//" data points"
                endif
                call this%ndffile%open(this%ndfname, "ua")
             endif
             
             if(coop_file_exists(this%chainname))then
                this%bestlike = coop_logZero
                this%bestparams = this%params
                call this%chainfile%open(this%chainname, "r")
                do while(this%chainfile%read_string(line))
                   read(line, *) this%knot
                   this%sum_mult = this%sum_mult + this%knot(1)
                   call this%chain%push(this%knot)
                   if(this%knot(2) .lt. this%bestlike)then
                      this%bestparams = this%knot(3:this%n+2)
                      this%bestlike = this%knot(2)
                   endif
                enddo
                call this%chainfile%close()
                if(this%chain%n.gt.0)then
                   this%params = this%knot(3:this%n+2)
                   this%fullparams(this%used) = this%params
                   this%loglike = this%knot(2)
                   this%mult = 1
                   this%loglike_is_exact = .false.
                   call this%chainfile%open(this%chainname, "a")
                   write(*,*) "continuing "//COOP_STR_OF(this%chain%n)//" lines from file "//trim(this%chainname)
                   if(.not. this%do_memsave) call this%chain%init()
                   if(associated(this%cosmology))then
                      call this%set_cosmology()
                   endif
                   this%derived_params = this%derived()
                else
                   this%do_overwrite = .true.
                endif
             else
                this%do_overwrite = .true.
             endif
          endif
          if(this%do_overwrite)then
             this%time = coop_systime_sec(0.d0)                          
             call this%chainfile%open(this%chainname, "w")
             if(associated(this%cosmology))then
                call this%set_cosmology()
             endif
             this%derived_params = this%derived()
             this%mult = 1.d0
             this%loglike = pool%LogLike(this)
             this%num_exact_calc = this%num_exact_calc + 1
             this%loglike_is_exact = .true.
             this%bestparams = this%params
             this%bestlike = this%loglike
          else
             this%time = coop_systime_sec(this%update_seconds*min(0.9d0, this%chain%n/200.d0))
          endif
       endif
       this%params_saved = this%params
       this%cosmology_saved => this%cosmology
       cosmology_changed = .false.
       if(this%do_fastslow)then
          if(this%loglike_is_exact)then
             if(this%cycl%next() .gt. this%n_slow)then
                vec(1:this%n_fast) = this%propose_fast_vec()*this%proposal_r()
                this%params(this%index_fast_start:this%n) = this%params_saved(this%index_fast_start:this%n) + matmul(this%mapping_fast, vec(1:this%n_fast))
                this%slow_changed = .false.
             else
                vec(1:this%n_slow) = this%propose_slow_vec()*this%proposal_r()
                this%params = this%params_saved + matmul(this%mapping(:, 1:this%n_slow), vec(1:this%n_slow))
                this%slow_changed = .true.             
             endif
          else
             vec(1:this%n_slow) = this%propose_slow_vec()*this%proposal_r()
             this%params = this%params_saved + matmul(this%mapping(:, 1:this%n_slow), vec(1:this%n_slow))
             this%slow_changed = .true.             
          endif            
       else
          vec = this%propose_vec()*this%proposal_r()
          this%params = this%params_saved + matmul(this%mapping, vec)
          this%slow_changed = .true.          
       endif
       this%fullparams(this%used) = this%params
       
       if(this%priorlike() .lt. coop_logZero)then
          if(this%do_ndf .and. this%slow_changed .and. this%like_approx%n .gt. this%n * 50 .and. this%converge_R .lt. 30.d0  .and. coop_random_unit().lt. this%approx_frac .and. (.not. this%is_drift) )then
             this%loglike_proposed = this%like_approx%eval(this%params)
             this%num_approx_calc = this%num_approx_calc + 1
             this%loglike_proposed_is_exact = .false.
          else
             if(associated(this%cosmology) .and. this%slow_changed)then
                this%cosmology => cosmo
                cosmology_changed = .true.
                call this%set_cosmology()
             endif
             this%loglike_proposed = pool%loglike(this)
             this%num_exact_calc = this%num_exact_calc + 1
             this%loglike_proposed_is_exact = .true.             
          endif
       else
          this%loglike_proposed = coop_logZero
       endif
       
       if(this%loglike_proposed .lt. coop_logZero .and. ((this%loglike_proposed - this%loglike)/this%temperature .lt. coop_random_exp() .or. this%is_drift))then
          this%accept = this%accept + 1
          this%knot(1) = this%mult
          this%knot(2) = this%loglike
          this%knot(3:this%n+2) = this%params_saved
          if(cosmology_changed)then
             this%cosmology => this%cosmology_saved
             this%cosmology = cosmo
          endif
          this%sum_mult = this%sum_mult + this%knot(1)
          if(this%do_memsave) &          
               call this%chain%push(this%knot)          
          if(this%chainfile%unit .ne. 0)then
             if(this%do_write_chain)then
                write(this%chainfile%unit, trim(this%form)) this%knot, this%derived_params
                if(this%do_flush)call flush(this%chainfile%unit)
             endif
          endif
          if(this%ndffile%unit .ne. 0 .and. this%loglike_is_exact .and. this%do_ndf)then
             write(this%ndffile%unit) this%knot
             if(this%do_flush)call flush(this%ndffile%unit)                
          endif
          this%loglike = this%loglike_proposed
          this%loglike_is_exact = this%loglike_proposed_is_exact
          
          this%derived_params = this%derived()          
          this%mult  = 1.d0
          if(this%loglike .lt. this%bestlike)then
             this%bestlike = this%loglike
             this%bestparams = this%params
          endif
       else
          if(this%ndffile%unit .ne. 0 .and. this%do_ndf .and. this%loglike_proposed_is_exact )then
             this%knot(1) = 0.01d0
             this%knot(2) = this%loglike_proposed
             this%knot(3:this%n+2) = this%params
             write(this%ndffile%unit) this%knot
             if(this%do_flush)call flush(this%ndffile%unit)
          endif
          if(cosmology_changed)then
             this%cosmology => this%cosmology_saved
#if DO_EFT_DE
             if(this%cosmology%f_alpha_M%initialized) &
                  call coop_eft_DE_set_Mpsq(this%cosmology%f_alpha_M)
#endif             
          endif
          this%params = this%params_saved
          this%fullparams(this%used) = this%params       
          this%reject = this%reject + 1
          this%mult = this%mult + 1.d0
       endif
    class default
       stop "for compatibility with old versions of gfortran McMC_step only support type coop_MCMC_params"
    end select

    if(.not. this%do_general_loglike .and. this%feedback .gt. 0)then
       if(mod(this%accept+this%reject, 9/this%feedback+1).eq. 0)then
          write(*,*) "on Node "//COOP_STR_OF(this%proc_id)//": step "//COOP_STR_OF(this%accept + this%reject)//", likelihood = "//COOP_STR_OF(this%loglike)//", accept ratio = "//COOP_STR_OF(dble(this%accept)/max(this%accept+this%reject,1)), " exact likelihoods ratio = "//COOP_STR_OF(dble(this%num_exact_calc)/max(this%num_exact_calc + this%num_approx_calc,1))

       endif
    endif
    this%step = this%step + 1        
    call this%update_propose()
  end subroutine coop_MCMC_params_MCMC_step

  subroutine coop_MCMC_params_findbest(this, Pool, temperature)
    class(coop_MCMC_params)::this
    type(coop_data_pool)::pool
    COOP_REAL,optional::temperature
    COOP_REAL::mintemp
    COOP_INT i, nsteps
    this%temperature = 10.d0
    if(present(temperature))then
       mintemp = max(min(temperature, 1.d0), 1.d-10)
    else
       mintemp = 0.001d0
    endif
    if(this%chainfile%unit .ne. 0) call this%chainfile%close()
    if(this%ndffile%unit .ne. 0) call this%ndffile%close()
    do while(this%temperature .gt. 0.0001d0)
       nsteps = 50 + ceiling(this%temperature * 45)       
       do i=1, nsteps
          call this%mcmc_step(pool)
       enddo
       call this%update_propose()
       this%temperature = this%temperature*0.8d0       
       this%params = this%bestparams
       this%fullparams(this%used) = this%bestparams
       call this%chain%init()
    enddo
    this%temperature = 1.d0
  end subroutine coop_MCMC_params_findbest




  function coop_dataset_loglike(this, mcmc) result(loglike)
    class(coop_dataset)::this
    type(coop_MCMC_params)::mcmc
    COOP_REAL::loglike
    if(this%off)then
       loglike = 0.d0
       return
    endif
    loglike = coop_LogZero
  end function coop_dataset_loglike
  

  function coop_dataset_SN_Simple_loglike(this, mcmc) result(loglike)
    class(coop_dataset_SN_Simple)::this
    type(coop_MCMC_params)::mcmc
    COOP_REAL::loglike, h0mpc
    COOP_REAL::mu_theory(this%n), Mbar
    COOP_INT i
    if(this%n .eq. 0 .or. this%off)then
       logLike = 0.d0
       return
    endif
    if(associated(mcmc%cosmology))then
       h0mpc = mcmc%cosmology%H0Mpc()
       !$omp parallel do
       do i=1, this%n
          mu_theory(i) = 5.d0*log10(mcmc%cosmology%dL_of_z(this%z(i))/h0mpc)
       enddo
       !$omp end parallel do
    else
       h0mpc = 0.7/3000.d0  !!this does not change the final loglike
       mu_theory(1) = coop_integrate(drz, 0.d0, this%z(1))
       do i=2, this%n
          mu_theory(i) = mu_theory(i-1) + coop_integrate(drz, this%z(i-1), this%z(i))
       enddo
       !$omp parallel do
       do i=1, this%n
          mu_theory(i) = 5.d0*log10(coop_r_of_chi(mu_theory(i), MCMC_OMEGA_K) *(1.d0+this%z(i))/h0mpc)
       enddo
       !$omp end parallel do
    endif
    Mbar = sum((mu_theory - this%mu)*this%invdmusq)/this%suminvdmusq
    loglike = sum((mu_theory-this%mu-Mbar)*(mu_theory-this%mu+Mbar)*this%invdmusq)/2.d0
  contains

    function drz(z)
      COOP_REAL z, drz
      drz = 1.d0/sqrt(MCMC_OMEGA_M*(1.d0+z)**3 + MCMC_OMEGA_K*(1.d0+z)**2 + MCMC_OMEGA_LAMBDA*(1.d0+z)**(3.d0*(1.d0+MCMC_W + MCMC_WA))*exp(-3.d0*MCMC_WA*z/(1.d0+z)) )
    end function drz

  end function coop_dataset_SN_Simple_loglike


  function coop_dataset_CMB_simple_loglike(this, mcmc) result(loglike)
    class(coop_dataset_CMB_simple)::this
    type(coop_mcmc_params)::mcmc
    COOP_REAL::loglike, vec(3)
    if(associated(mcmc%cosmology))then
       if(.not. this%has_invcov)then
          this%invcov(1,1) = this%ombh2_sigma**2
          this%invcov(2,2) = this%omch2_sigma**2
          this%invcov(3,3) = this%theta_sigma**2
          this%invcov(1,2) = this%ombh2_sigma*this%omch2_sigma*this%corr_ombh2_omch2
          this%invcov(2,1) = this%invcov(1,2)
          this%invcov(1,3) = this%ombh2_sigma*this%theta_sigma*this%corr_ombh2_theta
          this%invcov(3,1) = this%invcov(1,3)
          this%invcov(2,3) = this%omch2_sigma*this%theta_sigma*this%corr_omch2_theta
          this%invcov(3,2) = this%invcov(2,3)
          call coop_sympos_inverse(3, 3, this%invcov)
          this%has_invcov = .true.
       endif
       vec = (/ mcmc%cosmology%ombh2 -this%ombh2_center, &
            mcmc%cosmology%omch2 - this%omch2_center,  &
            100.d0*mcmc%cosmology%cosmomc_theta() - this%theta_center /)
       loglike = dot_product(vec, matmul(this%invcov, vec))/2.d0
    else
       stop "cannot use compressed CMB likelihood for models without cosmology"
    endif
  end function coop_dataset_CMB_simple_loglike

  function coop_dataset_WL_LogLike(this, mcmc) result(loglike)
    class(coop_dataset_WL)::this
    type(coop_mcmc_params)::mcmc
    COOP_REAL::loglike
    COOP_INT :: i
    loglike = 0.d0
    if(.not. associated(this%WLlike))return
    if(.not.associated(mcmc%cosmology)) stop "WL like requires cosmology input"
    do i=1, size(this%WLlike)
       loglike = loglike + this%WLLike(i)%loglike(mcmc%cosmology)
       if(loglike .ge. coop_logZero) return
    enddo
  end function coop_dataset_WL_LogLike
  
  function coop_dataset_BAO_loglike(this, mcmc) result(loglike)
    class(coop_dataset_BAO)::this
    type(coop_mcmc_params)::mcmc
    COOP_REAL::loglike
    COOP_INT :: i
    loglike = 0.d0
    if(.not. associated(this%baolike))return
    if(.not.associated(mcmc%cosmology)) stop "BAO like requires cosmology input"
    do i=1, size(this%baolike)
       loglike = loglike + this%baolike(i)%loglike(mcmc%cosmology)
       if(loglike .ge. coop_logZero) return
    enddo
  end function coop_dataset_BAO_loglike

  function coop_dataset_CMB_LogLike(this, mcmc) result(loglike)
    class(coop_dataset_CMB)::this
    type(coop_mcmc_params)::mcmc
    COOP_REAL::loglike
    COOP_REAL, dimension(:),allocatable::pars
    COOP_INT::i, inuis, ind, l
    loglike = 0.d0    
    if(.not. associated(this%cliklike))return
    if(.not. associated(mcmc%cosmology))then
       stop "for CMB likelihood you need to initialize the cosmology object"
    endif
    inuis = 1
    do i = 1, size(this%cliklike)
       if(.not. this%cliklike(i)%initialized) cycle
       if(this%cliklike(i)%numnames .gt. 0)then
          do inuis = 1, this%cliklike(i)%numnames
             ind = mcmc%index_of(trim(this%cliklike(i)%names(inuis)))
             if(ind .eq. 0)then
                write(*,*) "param["//trim(this%cliklike(i)%names(inuis))//"] not found"
                stop
             endif
             this%cliklike(i)%pars(inuis) = mcmc%fullparams(ind)
          enddo
          call this%cliklike(i)%set_cl_and_pars(mcmc%Cls_lensed, this%cliklike(i)%pars)          
       else
          call this%cliklike(i)%set_cl_and_pars(mcmc%Cls_lensed)                
       endif
       LogLike = LogLike + this%cliklike(i)%LogLike()
    enddo
  end function coop_dataset_CMB_LogLike

  function coop_dataset_HST_logLike(this, mcmc) result(loglike)
    class(coop_dataset_HST)::this
    type(coop_mcmc_params)::Mcmc
    COOP_REAL::loglike
    if(associated(this%HSTlike))then
       if(associated(mcmc%cosmology))then
          LogLike = this%HSTLike%LogLike(mcmc%cosmology%dA_of_z(this%HSTLike%zeff)/mcmc%cosmology%H0Mpc())
       else
          stop "For HST likelihood you need to initialize cosmology"
       endif
    else
       LogLike = 0.d0       
    endif
  end function coop_dataset_HST_logLike

  function coop_dataset_SN_JLA_loglike(this, mcmc) result(loglike)
    class(coop_dataset_SN_JLA)::this
    type(coop_mcmc_params)::mcmc
    COOP_REAL::loglike
    COOP_REAL grid_best, zhel, zcmb, alpha, beta
    COOP_INT grid_i, i, ind_alpha, ind_beta, j
   ! type(coop_file)::fp
    if(associated(this%JLALike))then
       if(.not. associated(mcmc%cosmology))stop "for JLA you need initialize cosmology"
       loglike = coop_logzero
       if(this%JLALike%n .eq. 0) stop "JLA data not loaded"
       !call fp%open("fidlum.txt", "r")
       do i=1, this%JLALike%n
          if(this%JLALike%sn(i)%has_absdist)then
             this%JLALike%lumdists(i) = 5.0*LOG10( this%JLALike%sn(i)%absdist)
          else
             zhel = this%JLALike%sn(i)%zhel
             zcmb = this%JLALike%sn(i)%zcmb
             this%JLALike%lumdists(i) = 5.d0*log10((1.0+zhel)/(1.0+zcmb) * mcmc%cosmology%dL_of_z(zcmb)/mcmc%cosmology%H0Mpc())
          endif
          !read(fp%unit,*) j, this%JLALike%lumdists(i)
          !if(j.ne.i) stop "Error in fidlum.txt"
       enddo
       !call fp%close()
       if (this%JLALike%marginalize) then
          !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC), PRIVATE(alpha,beta, grid_i)
          do grid_i = 1, this%JLALike%int_points
             alpha = this%JLALike%alpha_grid(grid_i)
             beta=this%JLALike%beta_grid(grid_i)
             this%JLALike%marge_grid(grid_i) = this%JLALike%alpha_beta_like(alpha, beta, this%JLALike%lumdists)
          end do

          grid_best = minval(this%JLALike%marge_grid,  mask=this%JLALike%marge_grid .lt. coop_logZero)
          loglike =  grid_best - log(sum(exp(-this%JLALike%marge_grid + grid_best),  &
               mask=this%JLALike%marge_grid .lt. coop_logZero) * this%JLALike%step_width_alpha*this%JLALike%step_width_beta)
       else
          ind_alpha = mcmc%index_of("alpha_JLA")
          ind_beta = mcmc%index_of("beta_JLA")
          if(ind_alpha .eq. 0)stop "param[alpha_JLA] is not found"
          if(ind_beta .eq. 0)stop "param[beta_JLA] is not found"          
          alpha = mcmc%fullparams(ind_alpha)
          beta = mcmc%fullparams(ind_beta)
          loglike =this%JLALike%alpha_beta_like(alpha, beta, this%JLALike%lumdists)
       end if
    else
       LogLike = 0.d0
    endif
  end function coop_dataset_SN_JLA_loglike

  function coop_Data_Pool_LogLike(this, mcmc) result(LogLike)
    class(coop_Data_Pool)this
    type(coop_mcmc_params)::mcmc
    COOP_INT::i
    COOP_REAL LogLike, tmp
    COOP_REAL,dimension(:,:),allocatable::Cls
    !!Prior
    LogLike = mcmc%PriorLike()    
    if(LogLike .ge. coop_LogZero) return
    if(associated(mcmc%cosmology))then
       if(mcmc%cosmology%h().lt. 0.01d0)then
          LogLike =coop_LogZero
          return
       endif
    endif    
    if(mcmc%do_general_loglike)then
       loglike = LogLike + mcmc%general_loglike()
       if(.not.(LogLike .lt. coop_LogZero)) return                          
    else
       !!Supernova
       if(associated(this%SN_Simple))then
          do i=1, size(this%SN_Simple)
             LogLike = LogLike + this%SN_Simple(i)%LogLike(mcmc)
             if(.not.(LogLike .lt. coop_LogZero)) return                   
          enddo
       endif

       !!simple CMB
       if(associated(this%CMB_Simple))then
          LogLike = LogLike + this%CMB_Simple%LogLike(mcmc)
          if(.not.(LogLike .lt. coop_LogZero)) return                     
       endif
       !!HST
       tmp = this%HST%LogLike(mcmc)
       if(mcmc%feedback .gt. 2 ) write(*,*) "HST like", tmp
       LogLike = LogLike + tmp
       if(.not.(LogLike .lt. coop_LogZero)) return
       !!JLA SN
       tmp = this%SN_JLA%LogLike(mcmc)
       LogLike = LogLike + tmp
       if(mcmc%feedback .gt. 2 ) write(*,*) "JLA like", tmp
       if(.not.(LogLike .lt. coop_LogZero)) return
       !!BAO
       tmp =  this%BAO%LogLike(mcmc)
       if(mcmc%feedback .gt. 2 ) write(*,*) "BAO like", tmp              
       LogLike = LogLike + tmp
       if(.not.(LogLike .lt. coop_LogZero)) return
       !!WL
       tmp = this%WL%LogLike(mcmc)
       if(mcmc%feedback .gt. 2) write(*,*) "WL Like", tmp
       LogLike = LogLike + tmp
       if(.not.(LogLike .lt. coop_LogZero)) return
       !!CMB
       tmp = this%CMB%LogLike(mcmc)
       if(mcmc%feedback .gt. 2 ) write(*,*) "CMB like", tmp                     
       LogLike = LogLike + tmp
       if(.not.(LogLike .lt. coop_LogZero)) return
    endif
  end function coop_Data_Pool_LogLike


  subroutine coop_dataset_SN_Simple_export(this, fname)
    class(coop_dataset_SN_Simple)::this
    COOP_UNKNOWN_STRING::fname
    type(coop_file)::fp
    COOP_INT i
    call fp%open(fname, "w")
    do i=1, this%n
       write(fp%unit, "(A, 3E16.7)") "SN ", this%z(i), this%mu(i), this%dmu(i)
    enddo
    call fp%close()

  end subroutine coop_dataset_SN_Simple_export

  subroutine coop_dataset_SN_Simple_import(this, fname)
    class(coop_dataset_SN_Simple)::this
    COOP_UNKNOWN_STRING::fname
    type(coop_file)::fp
    type(coop_list_realarr)::rl
    COOP_SINGLE p(3)
    COOP_INT i
    COOP_STRING line, name
    if(.not. coop_file_exists(fname))then
       write(*,*) trim(fname)//" does not exist"
       stop
    endif
    call rl%init()
    call fp%open(fname, "r")
    do
       read(fp%unit, "(A)", END=200, ERR=200) line
       line = adjustl(line)
       if(trim(line).eq."")cycle
       if(line(1:1).eq."#")cycle
       read(line, *) name, p   !!name, z, mu, dmu
       call rl%push(p)
    enddo
200 call fp%close()
    if(rl%n .eq. 0)then
       stop "dataset_SN_import: no data points"
    endif
    call rl%sort(1)
    if(allocated(this%z))deallocate(this%z, this%mu, this%dmu, this%invdmusq)
    this%n = rl%n
    allocate(this%z(rl%n), this%mu(rl%n), this%dmu(rl%n), this%invdmusq(rl%n))
    do i=1, rl%n
       p = rl%element(i)
       this%z(i) = p(1)
       this%mu(i) = p(2)
       this%dmu(i) = p(3)
       this%invdmusq(i) = 1.d0/(this%dmu(i)**2 + (this%pec_vel/this%z(i))**2)
    enddo
    this%suminvdmusq = sum(this%invdmusq)
    call rl%init()
  end subroutine coop_dataset_SN_Simple_import

  
  subroutine coop_dataset_SN_Simple_simulate(this, mcmc, z, dmu)
    class(coop_dataset_SN_Simple)::this
    type(coop_mcmc_params)::mcmc
    COOP_INT i, n
    COOP_REAL::h0mpc
    COOP_REAL,dimension(:),optional:: z, dmu
    if(present(z))then
       if(present(dmu))then
          n = size(z)
          if(n .ne. size(dmu))then
             stop "SN_simulate: z and dmu must of the same size"
          endif
          
          if(allocated(this%z))deallocate(this%z, this%mu, this%dmu, this%invdmusq)
          allocate(this%z(n), this%mu(n), this%dmu(n), this%invdmusq(n))
          this%z = z
          this%dmu = dmu
          this%invdmusq = 1.d0/(this%dmu**2 + (this%pec_vel/this%z)**2)
          this%suminvdmusq = sum(this%invdmusq)
       else
          stop "SN_simulate:  z and dmu must be passed simultaneously"
       endif
    endif

    if(associated(mcmc%cosmology))then
       h0mpc=mcmc%cosmology%H0Mpc()           
       !$omp parallel do
       do i=1, this%n
          this%mu(i) = 5.d0*log10(mcmc%cosmology%dL_of_z(z(i))/h0mpc)
       enddo
       !$omp end parallel do
    else
       h0mpc = 0.7d0/3000.d0
       this%mu(1) = coop_integrate(drz, 0.d0, this%z(1))
       do i=2, this%n
          this%mu(i) = this%mu(i-1) + coop_integrate(drz, this%z(i-1), this%z(i))
       enddo
       do i=1, this%n
          this%mu(i) = 5.d0*log10((1.d0+this%z(i))*coop_r_of_chi(this%mu(i), MCMC_OMEGA_K)/h0mpc)
       enddo
    endif
  contains

    function drz(z)
      COOP_REAL z, drz
      drz = 1.d0/sqrt(MCMC_OMEGA_M*(1.d0+z)**3 + MCMC_OMEGA_K*(1.d0+z)**2 + MCMC_OMEGA_LAMBDA*(1.d0+z)**(3.d0*(1.d0+MCMC_W + MCMC_WA))*exp(-3.d0*MCMC_WA*z/(1.d0+z)) )
    end function drz
  end subroutine coop_dataset_SN_Simple_simulate

  subroutine coop_MCMC_params_init(this,  ini)
    class(coop_MCMC_params)::this
    COOP_UNKNOWN_STRING::ini
    type(coop_file)::fp
    logical success
    COOP_REAL,dimension(:),allocatable::center, lower, upper, width, iniwidth, prior_sigma, prior_center
    COOP_INT i, iused, j, stat
    COOP_STRING::paramnames, prefix
    COOP_STRING val, lcname1, lcname2
    COOP_LONG_STRING::line
    COOP_INT::ncov
    type(coop_list_string)::sl
    COOP_REAL,dimension(:,:),allocatable::cov_read
    COOP_INT, dimension(:),allocatable::ind_read
    if(coop_file_exists(trim(ini)))then
       call coop_load_dictionary(ini, this%settings)
    else
       write(*, "(A)") "MCMC_params_init: cannot find ini file "//trim(ini)
       stop
    endif
  call coop_dictionary_lookup(this%settings, "action", this%action)
  if(trim(this%action) .eq. "") this%action = "TEST"  
    
    prefix = this%settings%value("chain_name")
    if(trim(prefix).eq."")then
       write(*,*) "MCMC_params_init: you need to have an entry chain_name = ... in the ini file "//trim(ini)
       stop
    endif
    call this%chain%init()
    call this%paramnames%free()
    call this%covmat%free()
    call this%cycl%free()
    call this%like_approx%free()
    this%prefix = trim(adjustl(prefix))    
    this%proc_id = coop_MPI_rank()
    this%num_exact_calc = 0
    this%num_approx_calc = 0
    this%accept = 0
    this%reject = 0
    this%step = 0
    call coop_dictionary_lookup(this%settings, "feedback", this%feedback, 1)
    call coop_dictionary_lookup(this%settings, "general_loglike", this%do_general_loglike, .false.)
    
    select case(trim(this%action))
    case("MCMC")
       call coop_dictionary_lookup(this%settings, "overwrite", this%do_overwrite, .false.)
       if(this%do_overwrite)then
          write(*,*) "Warning: overwriting chains when overwrite options is on"
          this%converge_R = coop_logZero       
       else
          if(coop_file_exists(trim(this%prefix)//".converge_stat"))then
             call fp%open(trim(this%prefix)//".converge_stat", "r")
             read(fp%unit, iostat = stat)  this%converge_R
             call fp%close()
             if(stat .ne. 0)this%converge_R = coop_logZero 
          else
             this%converge_R = coop_logZero          
          endif
       endif
       call coop_dictionary_lookup(this%settings, "update_seconds", this%update_seconds, 0)
       call coop_dictionary_lookup(this%settings, "update_steps", this%update_steps, 0)
       if(this%update_seconds .eq. 0 .and. this%update_steps .eq. 0)then
          this%update_seconds = 100000000 !!just do not update
       endif
       if(this%update_seconds .gt. 3600*24*30)then
          this%do_memsave = .false.
       else
          this%do_memsave = .true.
       endif
       call coop_dictionary_lookup(this%settings, "total_steps", this%total_steps, 50000)    
       call coop_dictionary_lookup(this%settings, "approx_frac", this%approx_frac, 0.d0)
       this%do_ndf  = (this%approx_frac .gt. 0.d0)
       if(this%do_ndf .and. this%feedback .ge.1  .and. this%proc_id.eq.0)write(*,"(A, G15.4)") "approximation fraction:", this%approx_frac
       call coop_dictionary_lookup(this%settings, "drift_frac", this%drift_frac, 0.d0)
       call coop_dictionary_lookup(this%settings, "drift_step", this%drift_step, 0.01d0)    
       this%do_drift = (this%drift_frac .gt. 0.d0)
       if(this%do_drift .and. this%feedback .ge. 1 .and. this%proc_id.eq.0)then
          write(*, "(A, G15.4)") "random drift probability: ", this%drift_frac
       endif
    end select
    if(this%do_general_loglike)then
       if(associated(this%cosmology))nullify(this%cosmology)
    endif
    if(associated(this%cosmology))then
       this%n_derived = coop_n_derived_with_cosmology
    else
       if(this%do_general_loglike)then
          this%n_derived = 0
       else
          this%n_derived = coop_n_derived_without_cosmology
       endif
    endif
    if(allocated(this%derived_params))deallocate(this%derived_params)
    allocate(this%derived_params(this%n_derived))
    paramnames = this%settings%value("paramnames")
    if(trim(paramnames).eq."")then
       write(*,*) "MCMC_params_init: you need to have an entry paramnames = ... in the ini file "//trim(ini)
       stop
    endif
    if(coop_file_exists(trim(paramnames)))then
       this%fulln  = coop_file_numlines(paramnames)
       call coop_load_dictionary(paramnames, this%paramnames, col_key = 1)
    else
       write(*, "(A)") "MCMC_params_init: cannot find file "//trim(paramnames)
       stop
    endif
    if(allocated(this%fullparams))deallocate(this%used, this%map2used, this%fullparams, this%params, this%lower, this%upper, this%center, this%width, this%mapping, this%params_saved, this%name, this%tex, this%knot, this%bestparams, this%prior_sigma, this%prior_center, this%has_prior, this%vecs)
    this%index_propose = 1
    allocate(this%map2used(this%fulln), this%fullparams(this%fulln), this%name(this%fulln), this%tex(this%fulln), lower(this%fulln), upper(this%fulln), width(this%fulln), iniwidth(this%fulln), center(this%fulln), prior_sigma(this%fulln), prior_center(this%fulln))
    do i= 1, this%fulln
       this%name(i) = trim(this%paramnames%key(i))
       this%tex(i) = trim(this%paramnames%val(i))
       call coop_dictionary_lookup(this%settings, "prior["//trim(this%name(i))//"]", val)
       if(trim(val).eq."")then
          prior_center(i) = 0.d0
          prior_sigma(i) = 0.d0
       else
          read(val, *) prior_center(i), prior_sigma(i)
       endif
       call coop_dictionary_lookup(this%settings, "param["//trim(this%name(i))//"]", val)
       if(trim(val).eq."")then
          write(*,*) "key  param["//trim(this%name(i))//"] is not found in the ini file "//trim(ini)
          stop
       elseif(scan(trim(adjustl(val)), " "//coop_tab) .eq. 0)then
          read(val, *) center(i)
          lower(i) = center(i)
          upper(i) = center(i)
          width(i) = 0.d0
          iniwidth(i) = 0.d0
       else
          read(val, *) center(i), lower(i), upper(i), width(i), iniwidth(i)
          if(lower(i) .gt. center(i) .or. upper(i) .lt.  center(i))then
             write(*,*) "mcmc_params_init: "//trim(this%name(i))//" boundary check failed"
             write(*,*) center(i), lower(i), upper(i), width(i), iniwidth(i)
             stop
          endif
          if(width(i) .eq. 0.d0)then
             lower(i) = center(i)
             upper(i) = center(i)
             iniwidth(i) = 0.d0
          elseif(upper(i) - lower(i) .eq. 0.d0)then
             lower(i) = center(i)
             upper(i) = center(i)
             width(i) = 0.d0
             iniwidth(i) = 0.d0
          endif
       endif
       if(iniwidth(i) .ne. 0.d0 .and. trim(this%action).ne. "TEST")then
          this%fullparams(i) = center(i) + iniwidth(i)*coop_random_Gaussian()
          do while(this%fullparams(i) .gt. upper(i) .or. this%fullparams(i) .lt. lower(i))
             this%fullparams(i) = center(i) + iniwidth(i)*coop_random_Gaussian()
          enddo
       else
          this%fullparams(i) = center(i)
       endif
    enddo
    this%n = count(width .ne. 0.d0)


    call coop_dictionary_lookup(this%settings,  "n_lc_priors", this%n_lc_priors, 0)
    if(this%n_lc_priors .gt. 0)then
       if(allocated(this%lc_prior_center))then
          deallocate(this%lc_prior_center, this%lc_prior_sigma, this%lc_prior_weight, this%lc_prior_index)
       endif
       allocate(this%lc_prior_center(this%n_lc_priors), this%lc_prior_sigma(this%n_lc_priors), this%lc_prior_weight(2,this%n_lc_priors), this%lc_prior_index(2,this%n_lc_priors))
       do i=1, this%n_lc_priors
          call coop_dictionary_lookup(this%settings, "lc_prior"//COOP_STR_OF(i), val)
          read(val, *) lcname1, lcname2, this%lc_prior_weight(:, i), this%lc_prior_center(i), this%lc_prior_sigma(i)
          this%lc_prior_index(1, i) = this%index_of(trim(adjustl(lcname1)))
          if( this%lc_prior_index(1, i) .eq. 0)then
             write(*,"(A)") trim(adjustl(lcname1))
             stop "parameter not found"
          endif
          this%lc_prior_index(2, i) = this%index_of(trim(adjustl(lcname2)))
          if( this%lc_prior_index(2, i) .eq. 0)then
             write(*,"(A)") trim(adjustl(lcname2))
             stop "parameter not found"
          endif          
       enddo
    endif

    
    this%form = "("//COOP_STR_OF(this%n+2+this%n_derived)//"E16.7)"
    allocate(this%used(this%n), this%params(this%n), this%lower(this%n), this%upper(this%n), this%center(this%n), this%width(this%n), this%mapping(this%n, this%n), this%params_saved(this%n), this%knot(this%n+2), this%bestparams(this%n), this%prior_sigma(this%n), this%prior_center(this%n), this%has_prior(this%n), this%vecs(this%n, this%n))
    call this%covmat%alloc(this%n)
    call this%cycl%init(this%n)
    iused= 0
    this%map2used = 0
    do i = 1, this%fulln
       if(width(i) .ne. 0.d0)then
          iused = iused + 1
          this%used(iused) = i
          this%map2used(i) = iused
       endif
    enddo
    this%params = this%fullparams(this%used)
    this%bestparams = this%params
    this%center = center(this%used)
    this%width = width(this%used)
    this%upper = upper(this%used)
    this%lower = lower(this%used)
    this%prior_sigma = prior_sigma(this%used)
    this%prior_center = prior_center(this%used)
    this%has_prior = (this%prior_sigma .gt. 0.d0)
    this%any_prior = any(this%has_prior)
    if(this%any_prior .and. this%feedback .ge. 1)then
       do i=1, this%n
          if(this%has_prior(i))then
             write(*,*) "prior["//trim(this%name(this%used(i)))//"] = "//COOP_STR_OF(this%prior_center(i))//" +/- "//COOP_STR_OF(this%prior_sigma(i))
          endif
       enddo
    endif
    if(trim(this%action)  .ne. "TEST")then
       if((.not. this%do_overwrite) .and. coop_file_exists(trim(this%prefix)//".runcov"))then
          call this%covmat%import(trim(this%prefix)//".runcov")
       else
          call this%covmat%diagonal(this%width)
          call coop_dictionary_lookup(this%settings, "propose_matrix", val)
          if(trim(val).ne."")then
             call fp%open(val, "r")
             read(fp%unit, "(A)") line
             if(line(1:1).eq."#")line = adjustl(line(2:))
             call coop_string_to_list(line, sl)
             ncov = sl%n
             allocate(cov_read(ncov, ncov), ind_read(ncov))
             call coop_read_matrix(fp%unit, cov_read, ncov, ncov, success)
             call fp%close()
             if(success)then
                do i = 1, ncov
                   ind_read(i) = this%paramnames%index(trim(sl%element(i)))
                   if(ind_read(i).ne.0) ind_read(i) = this%map2used(ind_read(i))
                enddo
                if(any(ind_read .ne. 0))then
                   do i=1, ncov
                      if(ind_read(i).eq.0)cycle
                      do j=1, ncov
                         if(ind_read(j).eq.0)cycle
                         this%covmat%c(ind_read(i), ind_read(j)) = cov_read(i, j)/this%covmat%sigma(ind_read(i))/this%covmat%sigma(ind_read(j))
                      enddo
                   enddo
                   call this%covmat%normalize()             
                   if(this%proc_id.eq.0 .and. this%feedback .ge. 1)write(*, *) "propose matrix "//trim(val)//" is loaded with "//COOP_STR_OF(count(ind_read.ne.0))//" usable parameters from totally "//COOP_STR_OF(ncov)//" parameters"
                endif
             else
                if(this%proc_id.eq.0)write(*,*) "propose matrix "//trim(val)//" is broken"
                stop
             endif
             deallocate(cov_read, ind_read)
             call sl%init()
          endif
       endif
       do i=1, this%n
          this%mapping(i,:) = this%covmat%L(i,:)*this%covmat%sigma(i)
       enddo


       this%n_fast = 0    
       call coop_dictionary_lookup(this%settings, "last_slow", val)
       if(trim(val).ne."")then
          this%n_slow = this%paramnames%index(trim(val))
          if(this%n_slow .ne. 0)then
             if(this%n_slow .lt. this%fulln)then
                this%n_fast = count(width(this%n_slow+1:this%fulln).gt.0.d0)
             endif
          endif
       endif
       this%do_fastslow = (this%n_fast .gt. 0)
       this%n_slow = this%n - this%n_fast
       this%index_fast_start = this%n_slow + 1
       this%index_propose_fast = 1
       this%index_propose_slow = 1

       if(this%do_fastslow)then
          if(allocated(this%mapping_fast))deallocate(this%mapping_fast)
          if(allocated(this%vecs_fast))deallocate(this%vecs_fast)
          if(allocated(this%vecs_slow))deallocate(this%vecs_slow)       
          allocate(this%mapping_fast(this%n_fast, this%n_fast))
          allocate(this%vecs_fast(this%n_fast, this%n_fast))
          allocate(this%vecs_slow(this%n_slow, this%n_slow))       
          this%mapping_fast = this%mapping(this%index_fast_start:this%n, this%index_fast_start:this%n)
       endif
    endif
    if(this%proc_id .eq. 0)then
       
       call coop_export_dictionary(trim(this%prefix)//".inputparams", this%settings)
       if(trim(this%action) .ne. "TEST")then
          call fp%open(trim(this%prefix)//".ranges", "w")
          do i=1, this%fulln
             write(fp%unit, "(A16, 2E16.7)") this%name(i), lower(i), upper(i)
          enddo
          call fp%close()
          call fp%open(trim(this%prefix)//".paramnames", "w")
          do i=1, this%n
             write(fp%unit, "(A16, A32)")  this%name(this%used(i)), this%tex(this%used(i))
          enddo
          if(associated(this%cosmology))then
             write(fp%unit, "(A16, A32)") "H0              ", "H_0      "                    
             write(fp%unit, "(A16, A32)") "omegam          ", "\Omega_m  "          
             write(fp%unit, "(A16, A32)") "omegal          ", "\Omega_\Lambda  "
             write(fp%unit, "(A16, A32)") "sigma8          ", "\sigma_8"
          else
             write(fp%unit, "(A16, A32)") "omegal          ", "\Omega_\Lambda  "
          endif
          call fp%close()
       endif
    endif
    deallocate(center, lower, upper, width, iniwidth, prior_sigma, prior_center)    

    !!load all the indices
    this%index_ombm2h2 = this%index_of("ombm2h2")
    if(this%index_ombm2h2.eq.0) &
       this%index_ombm2h2 = this%index_of("ombh2")
    if(this%index_ombm2h2.eq.0) &
         this%index_ombm2h2 = this%index_of("omegabh2")
    if(this%index_ombm2h2 .eq. 0) write(*,*) "Warning: cannot find the default parameter omega_b M^2 h^2"
    this%index_omcm2h2 = this%index_of("omcm2h2")
    if(this%index_omcm2h2 .eq.0 ) &    
         this%index_omcm2h2 = this%index_of("omch2")
    if(this%index_omcm2h2 .eq.0 ) &
       this%index_omcm2h2 = this%index_of("omegach2")
    if(this%index_omcm2h2 .eq. 0) write(*,*) "Warning: cannot find the default parameter omega_c M^2 h^2"
    
    this%index_theta = this%index_of("theta")
    this%index_tau = this%index_of("tau")
    this%index_logAm2tau = this%index_of("logAm2tau")
    this%index_logA = this%index_of("logA")
    this%index_ns = this%index_of("ns")
    this%index_mnu = this%index_of("mnu")                        
    this%index_nrun = this%index_of("nrun")
    this%index_r = this%index_of("r")
    this%index_nt =      this%index_of("nt")

!!dark energy equation of state parameters    
    this%index_de_w = this%index_of("de_w")
    this%index_de_wa = this%index_of("de_wa")
    this%index_de_epss = this%index_of("de_epss")
    this%index_de_epsinf = this%index_of("de_epsinf")
    this%index_de_zetas = this%index_of("de_zetas")
    this%index_de_betas = this%index_of("de_betas")                

#if DO_EFT_DE    
    this%index_de_alpha_M0 = this%index_of("de_alpha_M0")
    this%index_de_alpha_H0 = this%index_of("de_alpha_H0")
    this%index_de_alpha_K0 = this%index_of("de_alpha_K0")
    this%index_de_alpha_T0 = this%index_of("de_alpha_T0")    
    this%index_de_alpha_B0 = this%index_of("de_alpha_B0")
    call coop_dictionary_lookup(this%settings, "w_is_background", this%w_is_background, .false.)
#elif DO_COUPLED_DE
    this%index_de_Q = this%index_of("de_Q")
#endif    
    this%index_h = this%index_of("h")
    this%index_omegam = this%index_of("omegam")
    this%index_omegab = this%index_of("omegab")    
    this%index_omegak = this%index_of("omegak")
  end subroutine coop_MCMC_params_init

  subroutine coop_MCMC_params_get_lmax_from_data(this, pool)
    class(coop_MCMC_params)::this
    type(coop_data_pool)::pool
    COOP_INT::i
    if(associated(pool%CMB%ClikLike))then  !!want cls
       this%lmax  = 0
       do i = 1, size(pool%CMB%ClikLike)
          if(pool%CMB%ClikLike(i)%is_lensing)then
             this%lmax = max(this%lmax, maxval(pool%CMB%ClikLike(i)%lensing_lmaxs))             
          else
             this%lmax = max(this%lmax, maxval(pool%CMB%ClikLike(i)%lmax))
          endif
       enddo
       this%lmax = this%lmax + 100 !!buffer for good lensing
       if(allocated(this%Cls_scalar))then
          deallocate(this%Cls_Scalar, this%Cls_tensor, this%Cls_lensed)
       endif
       allocate(this%Cls_scalar(coop_num_cls, 2:this%lmax), this%Cls_tensor(coop_num_cls, 2:this%lmax), this%Cls_lensed(coop_num_cls, 2:this%lmax))
    else
       this%lmax = 0
    endif
  end subroutine coop_MCMC_params_get_lmax_from_data

  function coop_MCMC_params_proposal_r(this) result(r)
    class(coop_MCMC_params)::this    
    COOP_REAL::r
    if(coop_random_unit() .lt. 1.d0/3.d0)then
       r = coop_random_exp()
    else
       r = sqrt((coop_random_Gaussian()**2 + coop_random_Gaussian()**2)/2.d0)
    endif
    r = this%proposal_length*r
    this%is_drift = .false.
    if(this%do_drift)then
       if(coop_random_unit() .lt. this%drift_frac)then
          this%is_drift = .true.
          r = r * this%drift_step
       endif
    endif
  end function coop_MCMC_params_proposal_r

  function coop_MCMC_params_propose_vec(this) result(vec)
    class(coop_MCMC_params)::this    
    COOP_REAL::vec(this%n)
    if(this%n .eq. 1)then
       vec = 1.d0
       return
    endif
    if(this%index_propose .le. 1 .or. this%index_propose .gt. this%n)then
       call coop_random_rotation(this%vecs, this%n)
       this%index_propose = 1
    endif
    vec = this%vecs(:, this%index_propose)
    if(this%index_propose .ge. this%n)then
       this%index_propose = 1
    else
       this%index_propose = this%index_propose + 1
    endif
  end function coop_MCMC_params_propose_vec

  function coop_MCMC_params_propose_fast_vec(this) result(vec)
    class(coop_MCMC_params)::this        
    COOP_REAL::vec(this%n_fast)
    if(this%n_fast .eq. 1)then
       vec = 1.d0
       return
    endif
    if(this%index_propose_fast .le. 1 .or. this%index_propose_fast .gt. this%n_fast)then
       call coop_random_rotation(this%vecs_fast, this%n_fast)
       this%index_propose_fast = 1
    endif
    vec = this%vecs_fast(:, this%index_propose_fast)
    if(this%index_propose_fast .ge. this%n_fast)then
       this%index_propose_fast = 1
    else
       this%index_propose_fast = this%index_propose_fast + 1
    endif
  end function coop_MCMC_params_propose_fast_vec


  function coop_MCMC_params_propose_slow_vec(this) result(vec)
    class(coop_MCMC_params)::this        
    COOP_REAL::vec(this%n_slow)
    if(this%n_slow .eq. 1)then
       vec = 1.d0
       return
    endif
    if(this%index_propose_slow .le. 1 .or. this%index_propose_slow .gt. this%n_slow)then
       call coop_random_rotation(this%vecs_slow, this%n_slow)
       this%index_propose_slow = 1
    endif
    vec = this%vecs_slow(:, this%index_propose_slow)
    if(this%index_propose_slow .ge. this%n_slow)then
       this%index_propose_slow = 1
    else
       this%index_propose_slow = this%index_propose_slow + 1
    endif
  end function coop_MCMC_params_propose_slow_vec
  

!!$++++++++++++++++  Define Your Own Likelihood ++++++++++++++++++++++++
!!$set "mcmc_general = T" in ini file
  function coop_MCMC_params_general_loglike(this) result(loglike)
    COOP_INT,parameter::nparams = 22
    class(coop_MCMC_params)::this
    COOP_INT:: i, j 
    COOP_REAL::loglike
    COOP_REAL,dimension(nparams, nparams), save::invcov
    COOP_REAL::p(25)
    type(coop_file)::fp
    logical,save::init = .true.
    COOP_STRING::line
    COOP_REAL,dimension(nparams)::fiducial = (/ 0.022, 0.11, 1.041, 0.08, 3.07, 0.96 , 123. , 56., 120., 2.3, 26., 8., 0.96, 0.31, 0.61, 1., 1., 0.32, 3.5, 0.27, 0.14, 3.06 /)
    if(init)then
       call fp%open("covmats/lcdm.covmat", "r")
       read(fp%unit, "(A)") line
       do i=1, nparams
          read(fp%unit, *, END=100, ERR=100) p
          invcov(1:nparams, i) = p(1:nparams)
       enddo
       call fp%close()
       call coop_sympos_inverse(nparams, nparams, invcov)
       init = .false.
    endif
    if(this%n .ne. nparams)then
       write(*,*) "general model # of paramemters: "//COOP_STR_OF(nparams)
       write(*,*) "ini file # of parameters: "//COOP_STR_OF(this%n)
       stop
    endif
    loglike =0.5d0*dot_product(matmul(invcov, this%params-fiducial), this%params-fiducial) + coop_random_Gaussian()/10.d0
    if(loglike .lt. -1.d-2) call coop_MPI_abort("error in covmat")
    return
100 stop "Error in file"    
  end function coop_MCMC_params_general_loglike


end module coop_forecast_mod  
