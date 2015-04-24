module coop_forecast_mod
  use coop_wrapper_firstorder
  use coop_clik_mod
  use coop_HSTlike_mod
  use coop_SNlike_JLA_mod
  implicit none
#include "constants.h"

#define MCMC_OMEGA_M mcmc%fullparams(1)
#define MCMC_OMEGA_K mcmc%fullparams(2)
#define MCMC_W    mcmc%fullparams(3)
#define MCMC_WA    mcmc%fullparams(4)
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
     COOP_REAL::zstar = 1089.d0
     COOP_REAL::R_center = 1.7488
     COOP_REAL::R_sigma = 0.0074
   contains
     procedure::set_R_center => coop_dataset_CMB_simple_set_R_center     
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


  type coop_Data_Pool
!!these are simulations     
     type(coop_dataset_SN_Simple),dimension(:), pointer::SN_Simple => null()
     type(coop_dataset_CMB_Simple), pointer::CMB_Simple => null()    
     type(coop_dataset_CMB)::CMB
     type(coop_dataset_HST)::HST
     type(coop_dataset_SN_JLA)::SN_JLA
     
   contains
     procedure::LogLike => coop_data_pool_LogLike
  end type coop_Data_Pool

  
  type coop_MCMC_params
     COOP_INT::proc_id = 0
     type(coop_cosmology_firstorder),pointer::cosmology => null()   
     type(coop_file)::chainfile
     type(coop_dictionary)::settings     
     COOP_STRING::prefix, chainname
     COOP_STRING::form
     logical::do_flush = .false.
     logical::do_write_chain = .true.
     logical::do_overwrite = .false.
     logical::do_fastslow = .false.
     logical::do_memsave = .false.
     COOP_INT::n_fast = 0
     COOP_INT::index_fast_start = 0
     COOP_INT::n_slow = 0
     COOP_INT::fast_steps = 0
     COOP_INT::fast_per_round = 5
     COOP_REAL::time = 0.d0
     COOP_INT:: n = 0
     COOP_INT:: fulln = 0
     COOP_INT:: n_derived = 0
     COOP_REAL:: proposal_length = 1.d0
     COOP_REAL::bestlike = coop_LogZero
     COOP_REAL::loglike = coop_LogZero
     COOP_REAL::loglike_proposed = coop_LogZero
     COOP_REAL::mult = 0.d0
     COOP_REAL::sum_mult = 0.d0
     COOP_REAL::temperature = 1.d0
     COOP_INT,dimension(:),allocatable::used
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
     COOP_REAL,dimension(:),allocatable::prior_sigma     
     COOP_REAL,dimension(:),allocatable::prior_center     
     COOP_REAL,dimension(:,:),allocatable::covmat
     COOP_REAL,dimension(:,:),allocatable::propose
     COOP_REAL,dimension(:),allocatable::params_saved
     COOP_SINGLE, dimension(:),allocatable::knot
     COOP_REAL,dimension(:,:),allocatable::propose_fast
     type(coop_list_realarr)::chain
     type(coop_dictionary)::paramnames
     COOP_INT::accept, reject
     !!index for parameters
     COOP_INT::index_ombh2 = 0
     COOP_INT::index_omch2 = 0          
     COOP_INT::index_theta = 0
     COOP_INT::index_tau = 0
     COOP_INT::index_mnu = 0
     COOP_INT::index_logA = 0
     COOP_INT::index_ns = 0
     COOP_INT::index_nrun = 0     
     COOP_INT::index_r = 0
     COOP_INT::index_nt = 0     
     COOP_INT::index_de_w = 0
     COOP_INT::index_de_wa = 0
     COOP_INT::index_de_Q = 0
     COOP_INT::index_de_tracking_n = 0
     COOP_INT::index_de_dUdphi = 0
     COOP_INT::index_de_dlnQdphi = 0
     COOP_INT::index_de_d2Udphi2 = 0
     COOP_INT::index_h = 0
   contains
     procedure::init => coop_MCMC_params_init
     procedure::MCMC_step => coop_MCMC_params_MCMC_step
     procedure::update_propose => coop_MCMC_params_update_Propose
     procedure::findbest => coop_MCMC_params_findbest
     procedure::derived => coop_MCMC_params_derived
     procedure::priorLike => coop_MCMC_params_priorLike
     procedure::Set_Cosmology => coop_MCMC_params_Set_Cosmology
     procedure::index_of => coop_MCMC_params_index_of
     
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
    Prior = coop_LogZero
    if(any(this%params .gt. this%upper) .or. any(this%params.lt.this%lower))return
    if(associated(this%cosmology))then
       if(this%cosmology%h().lt. 0.01d0)return
    endif
        
    Prior = sum(((this%params - this%prior_center)/this%prior_sigma)**2, mask = this%has_prior)/2.d0
  end function coop_MCMC_params_priorLike

  !!set up %cosmology from %fullparams
  subroutine coop_MCMC_params_Set_Cosmology(this)
    class(coop_MCMC_params)::this
    COOP_REAL, parameter::h_b_i = 0.4d0
    COOP_REAL, parameter::h_t_i =  1.d0
    COOP_REAL::h_t, h_b, h_m, theta_t, theta_b, theta_m, theta_want
    if(.not. associated(this%cosmology)) stop "MCMC_params_set_cosmology: cosmology not allocated"
    if(this%index_theta .ne. 0)then
       theta_want = this%fullparams(this%index_theta)
       h_t = h_t_i
       h_b = h_b_i
       call calc_theta(h_t, theta_t)
       call calc_theta(h_b, theta_b)
       if(theta_t .lt. theta_want)then
          call this%cosmology%set_h(0.d0)
          return
       endif
       if(theta_b .gt. theta_want)then
          call this%cosmology%set_h(0.d0)
          return
       endif
       do while(h_t - h_b .gt. 1.d-4)
          h_m = (h_t + h_b)/2.d0
          call calc_theta(h_m, theta_m)
          if(theta_m .gt. theta_want)then
             h_t = h_m
             theta_t = theta_m
          else
             h_b = h_m
             theta_b = theta_m
          endif
       enddo
       if(abs(theta_m - theta_want).gt.1.d-7)then
          call calc_theta((h_t*(theta_want - theta_b)+h_b*(theta_t - theta_want))/(theta_t - theta_b), theta_m)
       endif
    elseif(this%index_h .ne. 0)then
       call setForH(this%fullparams(this%index_h))
    else
       stop "you need to use either theta or h for MCMC runs"
    endif
    call this%cosmology%setup_background()
    if(this%index_tau .ne. 0)then !!
       this%cosmology%optre = this%fullparams(this%index_tau)
       call this%cosmology%set_xe()
       if(this%index_logA .ne. 0)then
          this%cosmology%As = exp(this%fullparams(this%index_logA))*1.d-10
       else
          this%cosmology%As = 2.2d-9
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
       call this%cosmology%compute_source(0)
       if(this%cosmology%has_tensor)then
          call this%cosmology%compute_source(2)
       endif
    endif
  contains

    subroutine calc_theta(h, theta)
      COOP_REAL::h, theta
      call setForH(h)
      theta = this%cosmology%cosmomc_theta()*100.d0
    end subroutine calc_theta

    subroutine setforH(h)
      COOP_REAL::h
      COOP_REAL::Q, dlnQdphi, dUdphi, d2Udphi2, tracking_n, w, wa
      call this%cosmology%free()
      call this%cosmology%init(name = "Cosmology", id = 0, h = h)
      this%cosmology%ombh2 =this%fullparams(this%index_ombh2)
      this%cosmology%omch2 =this%fullparams(this%index_omch2)      
      !!baryon
      call this%cosmology%add_species(coop_baryon(this%cosmology%ombh2/h**2))
      !!radiation
      call this%cosmology%add_species(coop_radiation(this%cosmology%Omega_radiation()))
      !!neutrinos
      if(this%index_mnu .ne. 0)then
         if(this%fullparams(this%index_mnu).gt. 0.01d0)then  !!do massive neutrinos
            call this%cosmology%add_species(coop_neutrinos_massive( &
                 this%cosmology%Omega_nu_per_species_from_mnu_eV( this%fullparams(this%index_mnu) ) ,&
                 this%cosmology%Omega_massless_neutrinos_per_species()))
            call this%cosmology%add_species( coop_neutrinos_massless(this%cosmology%Omega_massless_neutrinos_per_species()*(this%cosmology%NNu()-1)))
         else
            call this%cosmology%add_species( coop_neutrinos_massless(this%cosmology%Omega_massless_neutrinos()))
         endif
      else
         call this%cosmology%add_species( coop_neutrinos_massless(this%cosmology%Omega_massless_neutrinos()))  
      endif
      if(this%index_de_tracking_n .ne. 0 .or. this%index_de_dUdphi .ne. 0 .or. this%index_de_Q .ne. 0)then  !!
         if(this%index_de_tracking_n .ne. 0)then
            tracking_n = this%fullparams(this%index_de_tracking_n)
         else
            tracking_n = 0.d0
         endif
         if(this%index_de_Q .ne. 0)then  !!coupled DE
            Q = this%fullparams(this%index_de_Q)
         else
            Q = 0.d0
         endif
         if(this%index_de_dlnQdphi .ne. 0)then
            dlnQdphi = this%fullparams(this%index_de_dlnQdphi)
         else
            dlnQdphi = 0.d0
         endif
         if(this%index_de_dUdphi .ne. 0)then
            dUdphi = this%fullparams(this%index_de_dUdphi)
         else
            dUdphi = 0.d0
         endif
         if(this%index_de_d2Udphi2 .ne. 0)then
            d2Udphi2 = this%fullparams(this%index_de_d2Udphi2)
         else
            d2Udphi2 = 0
         endif
         call coop_background_add_coupled_DE(this%cosmology, Omega_c = this%cosmology%omch2/h**2, Q = Q, tracking_n =  tracking_n, dlnQdphi = dlnQdphi, dUdphi = dUdphi, d2Udphi2 = d2Udphi2)         
      else
         call this%cosmology%add_species(coop_cdm(this%cosmology%omch2/h**2))
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
      endif
    end subroutine setforH

  end subroutine coop_MCMC_params_Set_Cosmology

  function  coop_MCMC_params_derived(this) result(derived)
    class(coop_MCMC_params)::this
    COOP_REAL:: derived(this%n_derived), phi, dlnVdphi, aeq, Q
    if(associated(this%cosmology))then
       derived(1) = this%cosmology%h()*100.d0
       derived(2) = this%cosmology%Omega_m
       derived(3) = 1.d0 - this%cosmology%Omega_m
       if(O0_DE(this%cosmology)%genre .eq. COOP_SPECIES_COUPLED)then
          aeq = 1.d0
          do while(O0_DE(this%cosmology)%density(aeq) .gt. O0_CDM(this%cosmology)%density(aeq))
             aeq = aeq*0.95
          enddo
          do while(O0_DE(this%cosmology)%density(aeq) .le. O0_CDM(this%cosmology)%density(aeq))
             aeq = aeq*1.01
          enddo
          do while(O0_DE(this%cosmology)%density(aeq) .gt. O0_CDM(this%cosmology)%density(aeq))
             aeq = aeq*0.995
          enddo
          phi = O0_DE(this%cosmology)%DE_phi(aeq)
          dlnVdphi = O0_DE(this%cosmology)%DE_dlnVdphi(phi)
          Q = O0_DE(this%cosmology)%DE_Q(aeq)
          derived(4) = (dlnVdphi + Q)**2/2.d0
       else
          derived(4) = 0.d0
       endif
    else
       derived(1) = 1.d0 - this%fullparams(1) - this%fullparams(2)
    endif
  end function coop_MCMC_params_derived

  subroutine coop_MCMC_params_update_Propose(this)
    class(coop_MCMC_params)::this
    COOP_INT::i, istart, i1, i2
    COOP_REAL::mean(this%n), cov(this%n, this%n), mult, diff(this%n)
    if(.not. this%do_memsave)stop "cannot update propose matrix when do_memsave is off"
    if(this%chain%n .lt. 20)return
    istart =  this%chain%n/4
    mult = 0.d0
    mean = 0.d0
    cov = 0.d0
    do i = istart, this%chain%n  !!only use 75% samples
       call this%chain%get_element(i, this%knot)
       mult  = mult + this%knot(1)
       mean = mean + this%knot(1)*this%knot(3:this%n+2)
    enddo
    mean = mean / mult
    do i = istart, this%chain%n  !!only use 75% samples
       call this%chain%get_element(i, this%knot)
       diff  = this%knot(3:this%n+2)  - mean
       do i1 = 1, this%n
          do i2 = 1, i1
             cov(i2, i1) = cov(i2, i1) + diff(i1)*diff(i2)*this%knot(1)
          enddo
       enddo
    enddo
    cov = cov/mult
    do i1 = 1, this%n
       do i2 = i1+1, this%n
          cov(i2, i1) = cov(i1, i2)
       enddo
    enddo
    this%covmat = cov
    this%propose = cov
    call coop_matsym_sqrt(this%propose)
    if(this%do_fastslow)then
       this%propose_fast = this%covmat(this%index_fast_start:this%n, this%index_fast_start:this%n)
       call coop_matsym_sqrt(this%propose_fast)
    endif
    
  end subroutine coop_MCMC_params_update_Propose


  

  subroutine coop_MCMC_params_MCMC_step(this, Pool)
    class(coop_MCMC_params)::this
    type(coop_data_pool)::pool
    COOP_REAL vec(this%n)
    COOP_INT numlines
    COOP_LONG_STRING::line
    if(this%n .le. 0) stop "MCMC: no varying parameters"
    select type(this)
    class is(coop_MCMC_params)
       if(this%accept+this%reject .eq. 0)then
          this%time = coop_systime_sec(.true.)
          call coop_random_init()
          call this%chain%init()
          this%sum_mult = 0.d0
          this%fast_steps = 0
          this%chainname = trim(this%prefix)//"_"//COOP_STR_OF(coop_MPI_Rank()+1)//".txt"
          if(.not. this%do_overwrite)then
             if(coop_file_exists(this%chainname))then
                this%bestlike = coop_logZero
                this%bestparams = this%params
                call this%chainfile%open(this%chainname, "r")
                do while(this%chainfile%read_string(line))
                   read(line, *) this%knot
                   this%sum_mult = this%sum_mult + this%knot(1)
                   if(this%do_memsave) &
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
                   call this%chainfile%open(this%chainname, "a")  
                else
                   this%do_overwrite = .true.
                endif
             else
                this%do_overwrite = .true.
             endif             
          endif
          if(this%do_overwrite)then
             call this%chainfile%open(this%chainname, "w")
             if(associated(this%cosmology))then
                call this%set_cosmology()
             endif
             this%loglike = pool%LogLike(this)
             this%bestparams = this%params
             this%bestlike = this%loglike
          endif
       endif
       this%params_saved = this%params
       if(this%do_fastslow .and. this%fast_steps .lt. this%fast_per_round)then
          this%fast_steps = this%fast_steps + 1
          call coop_random_get_Gaussian_vector(this%n_fast, vec(1:this%n_fast))
          vec(1:this%n_fast) =  vec(1:this%n_fast)*(this%proposal_length*coop_random_Gaussian())
          this%params(this%index_fast_start:this%n) = this%params_saved(this%index_fast_start:this%n) + matmul(this%propose_fast, vec(1:this%n_fast))          
       else
          this%fast_steps = 0
          call coop_random_get_Gaussian_vector(this%n, vec)
          vec =  vec*(this%proposal_length*coop_random_Gaussian())
          this%params = this%params_saved + matmul(this%propose, vec)          
       endif
       this%fullparams(this%used) = this%params
       if(associated(this%cosmology) .and. this%fast_steps.eq.0)then
          call this%set_cosmology()
       endif
       this%loglike_proposed = pool%loglike(this)
       if((this%loglike_proposed - this%loglike)/this%temperature .lt. coop_random_exp())then
          this%accept = this%accept + 1
          this%knot(1) = this%mult
          this%knot(2) = this%loglike
          this%knot(3:this%n+2) = this%params_saved
          this%sum_mult = this%sum_mult + this%knot(1)
          if(this%do_memsave) &          
               call this%chain%push(this%knot)          
          if(this%chainfile%unit .ne. 0)then
             if(this%do_write_chain)then
                write(this%chainfile%unit, trim(this%form)) this%knot, this%derived()
                if(this%do_flush .and. this%fast_steps.eq.0)call flush(this%chainfile%unit)
             endif
          endif
          this%loglike = this%loglike_proposed
          this%mult  = 1.d0
          if(this%loglike .lt. this%bestlike)then
             this%bestlike = this%loglike
             this%bestparams = this%params
          endif
       else
          this%params = this%params_saved
          this%fullparams(this%used) = this%params       
          this%reject = this%reject + 1
          this%mult = this%mult + 1.d0
       endif
    end select
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
    do while(this%temperature .gt. 0.0001d0)
       nsteps = 50 + ceiling(this%temperature * 45)       
       do i=1, nsteps
          call this%mcmc_step(pool)
       enddo
       call this%update_propose()
       this%params = this%bestparams
       this%fullparams(this%used) = this%bestparams
       this%temperature = this%temperature*0.8d0
       this%covmat = this%covmat*0.8d0
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
          mu_theory(i) = 5.d0*log10((1.d0+this%z(i))*mcmc%cosmology%luminosity_distance(1.d0/(1.d0+this%z(i)))/h0mpc)
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
    COOP_REAL::loglike
    if(associated(mcmc%cosmology))then
       loglike = ((mcmc%cosmology%distlss*sqrt(mcmc%cosmology%omega_m) - this%R_center)/this%R_sigma)**2/2.d0
    else
       loglike = ((coop_r_of_chi(coop_integrate(drz, 0.d0, this%zstar), MCMC_OMEGA_K)*sqrt(MCMC_OMEGA_M) - this%R_center)/this%R_sigma)**2/2.d0
    endif
    
  contains
    function drz(z)
      COOP_REAL z, drz
      drz = 1.d0/sqrt(MCMC_OMEGA_M*(1.d0+z)**3 + MCMC_OMEGA_K*(1.d0+z)**2 + MCMC_OMEGA_LAMBDA*(1.d0+z)**(3.d0*(1.d0+MCMC_W + MCMC_WA))*exp(-3.d0*MCMC_WA*z/(1.d0+z)) + 9.d-5*(1.d0+z)**4)  !!a fiducial Omega_r is put in here
    end function drz
    
  end function coop_dataset_CMB_simple_loglike


  subroutine coop_dataset_CMB_simple_set_R_center(this, Omega_m) 
    class(coop_dataset_CMB_simple)::this
    COOP_REAL Omega_m
    this%R_center = coop_integrate(drz, 0.d0, this%zstar)*sqrt(omega_m)
  contains
    function drz(z)
      COOP_REAL z, drz
      drz = 1.d0/sqrt(Omega_m*(1.d0+z)**3 + (1.d0-Omega_m) + 9.d-5*(1.d0+z)**4)  !!a fiducial Omega_r is put in here
    end function drz
    
  end subroutine coop_dataset_CMB_simple_set_R_center
  


  function coop_dataset_CMB_LogLike(this, mcmc) result(loglike)
    class(coop_dataset_CMB)::this
    type(coop_mcmc_params)::mcmc
    COOP_REAL::loglike
    COOP_REAL, dimension(:,:), allocatable::Cls_scalar, Cls_tensor, Cls_lensed
    COOP_REAL, dimension(:),allocatable::pars
    COOP_INT::lmax, i, inuis, ind, l
    loglike = 0.d0    
    if(.not. associated(this%cliklike))return
    if(.not. associated(mcmc%cosmology))then
       stop "for CMB likelihood you need to initialize the cosmology object"
    endif
    lmax = 0
    do i = 1, size(this%cliklike)
       lmax = max(lmax, maxval(this%cliklike(i)%lmax))
    enddo
    lmax = lmax + 100  !!use buffer to get better lensing cls
    allocate(Cls_Scalar(coop_num_Cls, 2:lmax), Cls_Tensor(coop_num_Cls, 2:lmax), Cls_lensed(coop_num_Cls, 2:lmax))
    call mcmc%cosmology%source(0)%get_all_cls(2, lmax, Cls_scalar)
    call coop_get_lensing_Cls(2, lmax, Cls_Scalar, Cls_lensed)
    Cls_lensed = Cls_lensed + Cls_scalar
    if(mcmc%cosmology%has_tensor)then
       call mcmc%cosmology%source(2)%get_all_cls( 2, lmax, Cls_tensor)
       Cls_lensed = Cls_lensed + Cls_tensor       
    endif
    Cls_lensed = Cls_lensed*((mcmc%cosmology%Tcmb())**2*1.d12)
    inuis = 1
    do i = 1, size(this%cliklike)
       if(this%cliklike(i)%numnames .gt. 0)then
          do inuis = 1, this%cliklike(i)%numnames
             ind = mcmc%index_of(trim(this%cliklike(i)%names(inuis)))
             if(ind .eq. 0)then
                write(*,*) "param["//trim(this%cliklike(i)%names(inuis))//"] not found"
                stop
             endif
             this%cliklike(i)%pars(inuis) = mcmc%fullparams(ind)
          enddo
          call this%cliklike(i)%set_cl_and_pars(Cls_lensed, this%cliklike(i)%pars)
       else
          call this%cliklike(i)%set_cl_and_pars(Cls_lensed)                
       endif
       LogLike = LogLike - this%cliklike(i)%LogLike()
    enddo
    deallocate(Cls_Scalar, Cls_tensor, Cls_lensed)
  end function coop_dataset_CMB_LogLike

  function coop_dataset_HST_logLike(this, mcmc) result(loglike)
    class(coop_dataset_HST)::this
    type(coop_mcmc_params)::Mcmc
    COOP_REAL::loglike
    if(associated(this%HSTlike))then
       if(associated(mcmc%cosmology))then
          LogLike = this%HSTLike%LogLike(mcmc%cosmology%angular_diameter_distance(1.d0/(1.d0+this%HSTLike%zeff))/mcmc%cosmology%H0Mpc())
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
             this%JLALike%lumdists(i) = 5.d0*log10((1.0+zhel)/(1.0+zcmb) * mcmc%cosmology%luminosity_distance(1.d0/(1.d0+zcmb))/mcmc%cosmology%H0Mpc())
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
    COOP_REAL LogLike
    COOP_REAL,dimension(:,:),allocatable::Cls
    !!Prior
    LogLike = mcmc%PriorLike()
    if(LogLike .ge. coop_LogZero) return
    !!Supernova
    if(associated(this%SN_Simple))then
       do i=1, size(this%SN_Simple)
          LogLike = LogLike + this%SN_Simple(i)%LogLike(mcmc)
          if(LogLike .ge. coop_LogZero) return          
       enddo
    endif
    !!simple CMB
    if(associated(this%CMB_Simple))then
       LogLike = LogLike + this%CMB_Simple%LogLike(mcmc)
       if(LogLike .ge. coop_LogZero) return                 
    endif
    LogLike = LogLike + this%SN_JLA%LogLike(mcmc)
    if(LogLike .ge. coop_LogZero) return    
    !!CMB
    LogLike = LogLike + this%CMB%LogLike(mcmc)
    if(LogLike .ge. coop_LogZero) return
    LogLike = LogLike + this%HST%LogLike(mcmc)
    if(LogLike .ge. coop_LogZero) return
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
          this%mu(i) = 5.d0*log10((1.d0+this%z(i))*mcmc%cosmology%luminosity_distance(1.d0/(1.d0+z(i)))/h0mpc)
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
    COOP_INT i, iused
    COOP_STRING::paramnames, prefix
    COOP_STRING val
    if(coop_file_exists(trim(ini)))then
       call coop_load_dictionary(ini, this%settings)
    else
       write(*, "(A)") "MCMC_params_init: cannot find ini file "//trim(ini)
       stop
    endif
    prefix = this%settings%value("chain_name")
    if(trim(prefix).eq."")then
       write(*,*) "MCMC_params_init: you need to have an entry chain_name = ... in the ini file "//trim(ini)
       stop
    endif
    call coop_dictionary_lookup(this%settings, "overwrite", this%do_overwrite, .false.)
    call this%chain%init()
    call this%paramnames%free()
    this%prefix = trim(adjustl(prefix))
    this%proc_id = coop_MPI_rank()
    if(associated(this%cosmology))then
       this%n_derived = coop_n_derived_with_cosmology
    else
       this%n_derived = coop_n_derived_without_cosmology       
    endif
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
    if(allocated(this%fullparams))deallocate(this%used, this%fullparams, this%params, this%lower, this%upper, this%center, this%width, this%covmat, this%propose, this%params_saved, this%name, this%tex, this%knot, this%bestparams, this%prior_sigma, this%prior_center, this%has_prior)
    allocate(this%fullparams(this%fulln), this%name(this%fulln), this%tex(this%fulln), lower(this%fulln), upper(this%fulln), width(this%fulln), iniwidth(this%fulln), center(this%fulln), prior_sigma(this%fulln), prior_center(this%fulln))
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
       if(iniwidth(i) .ne. 0.d0)then
          this%fullparams(i) = center(i) + iniwidth(i)*coop_random_Gaussian()
          do while(this%fullparams(i) .gt. upper(i) .or. this%fullparams(i) .lt. lower(i))
             this%fullparams(i) = center(i) + iniwidth(i)*coop_random_Gaussian()
          enddo
       else
          this%fullparams(i) = center(i)
       endif
    enddo
    this%n = count(width .ne. 0.d0)
    this%form = "("//COOP_STR_OF(this%n+2+this%n_derived)//"E16.7)"
    allocate(this%used(this%n), this%params(this%n), this%lower(this%n), this%upper(this%n), this%center(this%n), this%width(this%n), this%covmat(this%n, this%n), this%propose(this%n, this%n), this%params_saved(this%n), this%knot(this%n+2), this%bestparams(this%n), this%prior_sigma(this%n), this%prior_center(this%n), this%has_prior(this%n))
    iused= 0
    do i = 1, this%fulln
       if(width(i) .ne. 0.d0)then
          iused = iused + 1
          this%used(iused) = i
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
    call coop_dictionary_lookup(this%settings, "propose_matrix", val)
    if(trim(val).ne."")then
       call fp%open(val, "r")
       call coop_read_matrix(fp%unit, this%n, this%n, this%covmat, success)
       call fp%close()
       if(success)then
          this%propose = this%covmat
          call coop_matsym_sqrt(this%propose)
       else
          write(*,*) "propose matrix "//trim(val)//" is broken"
          stop
       endif
    elseif(coop_file_exists(trim(this%prefix)//".covmat"))then
       call fp%open(trim(this%prefix)//".covmat", "r")
       call coop_read_matrix(fp%unit, this%n, this%n, this%covmat, success)
       call fp%close()
       if(success)then
          this%propose = this%covmat
          call coop_matsym_sqrt(this%propose)
       endif
    else
       success = .false.       
    endif
    if(.not. success)then
       this%covmat = 0.d0
       this%propose = 0.d0
       do i=1, this%n
          this%covmat(i, i) = this%width(i)**2
          this%propose(i, i) = this%width(i)
       enddo
    endif
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
    if(this%do_fastslow)then
       this%n_slow = this%n - this%n_fast
       this%index_fast_start = this%n_slow + 1
       if(allocated(this%propose_fast))deallocate(this%propose_fast)
       allocate(this%propose_fast(this%n_fast, this%n_fast))
       this%propose_fast = this%covmat(this%n_slow+1:this%n, this%n_slow+1:this%n)
       call coop_matsym_sqrt(this%propose_fast)
    endif
    if(this%proc_id .eq. 0)then
       
       call coop_export_dictionary(trim(this%prefix)//".inputparams", this%settings)
       call fp%open(trim(this%prefix)//".ranges", "w")
       do i=1, this%fulln
          write(fp%unit, "(A32, 2E16.7)") this%name(i), lower(i), upper(i)
       enddo
       call fp%close()
       call fp%open(trim(this%prefix)//".paramnames", "w")
       do i=1, this%n
          write(fp%unit, "(2A16)")  this%name(this%used(i)), this%tex(this%used(i))
       enddo
       if(associated(this%cosmology))then
          write(fp%unit, "(2A16)") "H0              ", "H_0      "                    
          write(fp%unit, "(2A16)") "omegam          ", "\Omega_m  "          
          write(fp%unit, "(2A16)") "omegal          ", "\Omega_\Lambda  "
          write(fp%unit, "(2A16)") "epss            ", "\epsilon_s  "                    
       else
          write(fp%unit, "(2A16)") "omegal          ", "\Omega_\Lambda  "
       endif
       call fp%close()
       
    endif
    deallocate(center, lower, upper, width, iniwidth, prior_sigma, prior_center)    

    !!load all the indices
    this%index_ombh2 = this%index_of("ombh2")
    if(this%index_ombh2.eq.0)then
       this%index_ombh2 = this%index_of("omegabh2")
    endif
    this%index_omch2 = this%index_of("omch2")
    if(this%index_omch2 .eq.0 )then
       this%index_omch2 = this%index_of("omegach2")
    endif
    this%index_theta = this%index_of("theta")
    this%index_tau = this%index_of("tau")
    this%index_logA = this%index_of("logA")
    this%index_ns = this%index_of("ns")
    this%index_mnu = this%index_of("mnu")                        
    this%index_nrun = this%index_of("nrun")
    this%index_r = this%index_of("r")
    this%index_nt =      this%index_of("nt")
    this%index_de_w = this%index_of("de_w")
    this%index_de_wa = this%index_of("de_wa")
    this%index_de_Q = this%index_of("de_Q")
    this%index_de_tracking_n = this%index_of("de_tracking_n")
    this%index_de_dUdphi = this%index_of("de_dUdphi")
    this%index_de_dlnQdphi = this%index_of("de_dlnQdphi")
    this%index_de_d2Udphi2 = this%index_of("de_d2Udphi2")
    this%index_h = this%index_of("h")

    

  end subroutine coop_MCMC_params_init


end module coop_forecast_mod  
