module coop_forecast_mod
  use coop_wrapper_firstorder
#ifdef HAS_CLIK
  use clik
#endif  
  implicit none
#include "constants.h"

#define MCMC_OMEGA_M mcmc%fullparams(1)
#define MCMC_OMEGA_K mcmc%fullparams(2)
#define MCMC_W    mcmc%fullparams(3)
#define MCMC_WA    mcmc%fullparams(4)
#define MCMC_OMEGA_LAMBDA  (1.d0 - MCMC_OMEGA_M - MCMC_OMEGA_K)  

  
  COOP_REAL, parameter::coop_LogZero = 1.d30
  type coop_DataSet
     COOP_STRING::name
     type(coop_arguments)::args
     logical::off = .false.
     COOP_INT:: n_nuis = 0
     COOP_REAL, dimension(:), allocatable::nuis
     COOP_REAL, dimension(:), allocatable::nuis_upper
     COOP_REAL, dimension(:), allocatable::nuis_lower
     COOP_REAL, dimension(:), allocatable::nuis_center
     COOP_REAL, dimension(:), allocatable::nuis_sigma     
   contains
     procedure::prior => coop_dataset_prior     
     procedure::LogLike => coop_dataset_loglike
  end type coop_DataSet

  type, extends(coop_DataSet):: coop_DataSet_SN
     COOP_INT::n = 0
     COOP_REAL,dimension(:),allocatable::z, mu, dmu, invdmusq
     COOP_REAL::suminvdmusq
     COOP_REAL::pec_vel = 400.d0/3.d5
   contains
     procedure::LogLike =>coop_dataset_SN_loglike
     procedure::import => coop_dataset_SN_import
     procedure::simulate => coop_dataset_SN_simulate
     procedure::export => coop_dataset_SN_export   
  end type coop_DataSet_SN


  type coop_Data_Pool
     type(coop_DataSet_SN),dimension(:), pointer::SN
!     type(coop_DataSet_SN),dimension(:), pointer::BAO
!     type(coop_DataSet_SN),dimension(:), pointer::CMB
!     type(coop_DataSet_SN),dimension(:), pointer::WL
!     type(coop_DataSet_SN),dimension(:), pointer::Lya
!     type(coop_DataSet_SN),dimension(:), pointer::HST     
   contains
     procedure::Prior => coop_data_pool_prior
     procedure::LogLike => coop_data_pool_LogLike
  end type coop_Data_Pool

  
  type coop_MCMC_params
     COOP_INT::proc_id = 0
     type(coop_cosmology_firstorder),pointer::cosmology => null()   
     type(coop_file)::chainfile
     COOP_STRING::prefix
     COOP_STRING::form
     COOP_INT:: n = 0
     COOP_INT:: fulln = 0
     COOP_INT:: n_derived = 1
     COOP_REAL:: proposal_length = 1.2d0
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
     COOP_REAL,dimension(:,:),allocatable::covmat
     COOP_REAL,dimension(:,:),allocatable::propose
     COOP_REAL,dimension(:),allocatable::params_saved
     COOP_SINGLE, dimension(:),allocatable::knot
     type(coop_list_realarr)::chain
     COOP_INT::accept, reject
   contains
     procedure::init => coop_MCMC_params_init
     procedure::MCMC_step => coop_MCMC_params_MCMC_step
     procedure::update_propose => coop_MCMC_params_update_Propose
     procedure::findbest => coop_MCMC_params_findbest
     procedure::derived => coop_MCMC_params_derived
  end type coop_MCMC_params



contains

  function  coop_MCMC_params_derived(this) result(derived)
    class(coop_MCMC_params)::this
    COOP_REAL:: derived(this%n_derived)
    derived(1) = 1.d0 - this%fullparams(1) - this%fullparams(2)
  end function coop_MCMC_params_derived

  subroutine coop_MCMC_params_update_Propose(this)
    class(coop_MCMC_params)::this
    COOP_INT::i, istart, i1, i2
    COOP_REAL::mean(this%n), cov(this%n, this%n), mult, diff(this%n)
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
  end subroutine coop_MCMC_params_update_Propose


  subroutine coop_MCMC_params_init(this, prefix, paramnames, ini)
    class(coop_MCMC_params)::this
    COOP_UNKNOWN_STRING::prefix, ini, paramnames
    type(coop_dictionary)::dict, pn
    type(coop_file)::fp
    logical success
    COOP_REAL,dimension(:),allocatable::center, lower, upper, width, iniwidth
    COOP_INT i, iused
    COOP_STRING val
    this%prefix = trim(adjustl(prefix))
    this%proc_id = coop_MPI_rank()
    if(coop_file_exists(trim(paramnames)))then
       this%fulln  = coop_file_numlines(paramnames)
       call coop_load_dictionary(paramnames, pn, col_key = 1, col_value = 2)
    else
       write(*, "(A)") "MCMC_params_init: cannot find file "//trim(paramnames)
       stop
    endif
    if(coop_file_exists(trim(ini)))then
       call coop_load_dictionary(ini, dict)
    else
       write(*, "(A)") "MCMC_params_init: cannot find file "//trim(ini)
       stop
    endif
    if(allocated(this%fullparams))deallocate(this%used, this%fullparams, this%params, this%lower, this%upper, this%center, this%width, this%covmat, this%propose, this%params_saved, this%name, this%tex, this%knot, this%bestparams)
    allocate(this%fullparams(this%fulln), this%name(this%fulln), this%tex(this%fulln), lower(this%fulln), upper(this%fulln), width(this%fulln), iniwidth(this%fulln), center(this%fulln))
    do i= 1, this%fulln
       this%name(i) = trim(pn%key(i))
       this%tex(i) = trim(pn%val(i))
       call coop_dictionary_lookup(dict, "param["//trim(this%name(i))//"]", val)
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
    allocate(this%used(this%n), this%params(this%n), this%lower(this%n), this%upper(this%n), this%center(this%n), this%width(this%n), this%covmat(this%n, this%n), this%propose(this%n, this%n), this%params_saved(this%n), this%knot(this%n+2), this%bestparams(this%n))
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
    call coop_dictionary_lookup(dict, "propose_matrix", val)
    if(trim(val).ne."")then
       call fp%open(val, "r")
       call coop_read_matrix(fp%unit, this%n, this%n, this%propose, success)
       call fp%close()
       if(.not. success)then
          write(*,*) "propose matrix "//trim(val)//" is broken"
          stop
       endif
       this%covmat = matmul(this%propose, this%propose)
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
    if(this%proc_id .eq. 0)then
       
       call coop_export_dictionary(trim(this%prefix)//".inputparams", dict)
       call fp%open(trim(this%prefix)//".ranges", "w")
       do i=1, this%fulln
          write(fp%unit, "(A32, 2E16.7)") this%name(i), lower(i), upper(i)
       enddo
       call fp%close()
       call fp%open(trim(this%prefix)//".paramnames", "w")
       do i=1, this%n
          write(fp%unit, "(2A16)")  this%name(this%used(i)), this%tex(this%used(i))
       enddo
       write(fp%unit, "(2A16)") "omegal          ", "\Omega_\Lambda  "
       call fp%close()
       
    endif
    deallocate(center, lower, upper, width, iniwidth)    
    call dict%free()
    call pn%free()

  end subroutine coop_MCMC_params_init
  

  subroutine coop_MCMC_params_MCMC_step(this, Pool)
    class(coop_MCMC_params)::this
    type(coop_data_pool)::pool
    COOP_REAL vec(this%n)
    if(this%n .le. 0) stop "MCMC: no varying parameters"
    select type(this)
    class is(coop_MCMC_params)
       if(this%accept+this%reject .eq. 0)then
          call this%chainfile%open(trim(this%prefix)//"_"//COOP_STR_OF(this%proc_id+1)//".txt", "w")          
          call this%chain%init()
          this%loglike = pool%LogLike(this)
          this%bestparams = this%params
          this%bestlike = this%loglike
       endif
       this%params_saved = this%params
       vec =  coop_random_vector(this%n)*(this%proposal_length*coop_random_Gaussian())       
       this%params = this%params_saved + matmul(this%propose, vec)
       if(any(this%params .gt. this%upper) .or. any(this%params .lt. this%lower))then
          this%params = this%params_saved
          this%fullparams(this%used) = this%params       
          this%reject = this%reject + 1
          this%mult = this%mult + 1.d0
          return
       endif
       this%fullparams(this%used) = this%params           
       this%loglike_proposed = pool%loglike(this)
       if((this%loglike_proposed - this%loglike)/this%temperature .lt. coop_random_exp())then
          this%accept = this%accept + 1
          this%knot(1) = this%mult
          this%knot(2) = this%loglike
          this%knot(3:this%n+2) = this%params_saved
          call this%chain%push( this%knot)          
          if(this%chainfile%unit .ne. 0)then
             write(this%chainfile%unit, trim(this%form)) this%knot, this%derived()
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



  function coop_dataset_prior(this, mcmc) result(prior)
    class(coop_dataset)::this
    type(coop_MCMC_params)::mcmc
    COOP_REAL::prior
    if(this%n_nuis .eq. 0 .or. this%off)then
       prior = 0.d0
       return
    endif
    if(any(this%nuis .gt. this%nuis_upper) .or. any(this%nuis .lt. this%nuis_lower))then
       prior = coop_LogZero
       return
    endif
    prior = sum(((this%nuis - this%nuis_center)/this%nuis_sigma)**2)/2.d0
  end function coop_dataset_prior
  

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
  

  function coop_dataset_SN_loglike(this, mcmc) result(loglike)
    class(coop_dataset_SN)::this
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

  end function coop_dataset_SN_loglike

  function coop_Data_Pool_prior(this, mcmc) result(prior)
    class(coop_Data_Pool)this
    type(coop_mcmc_params)::mcmc
    COOP_REAL prior
    if(associated(mcmc%cosmology))then
       if(mcmc%cosmology%h() .le. 0.d0)then
          prior = coop_LogZero
          return
       endif
       prior = 0.d0
       return
    endif
    if(MCMC_OMEGA_LAMBDA .lt. 0.d0)then
       prior = coop_LogZero
    else
       prior = 0.d0
    endif
  end function coop_Data_Pool_prior

  function coop_Data_Pool_LogLike(this, mcmc) result(LogLike)
    class(coop_Data_Pool)this
    type(coop_mcmc_params)::mcmc
    COOP_INT::i
    COOP_REAL LogLike
    LogLike = this%Prior(mcmc)
    if(LogLike .ge. coop_LogZero) return
    do i=1, size(this%SN)
       LogLike = LogLike + this%SN(i)%LogLike(mcmc)
    enddo
  end function coop_Data_Pool_LogLike


  subroutine coop_dataset_SN_export(this, fname)
    class(coop_dataset_SN)::this
    COOP_UNKNOWN_STRING::fname
    type(coop_file)::fp
    COOP_INT i
    call fp%open(fname, "w")
    do i=1, this%n
       write(fp%unit, "(A, 3E16.7)") "SN ", this%z(i), this%mu(i), this%dmu(i)
    enddo
    call fp%close()

  end subroutine coop_dataset_SN_export

  subroutine coop_dataset_SN_import(this, fname)
    class(coop_dataset_SN)::this
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
  end subroutine coop_dataset_SN_import

  
  subroutine coop_dataset_SN_simulate(this, mcmc, z, dmu)
    class(coop_dataset_SN)::this
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
  end subroutine coop_dataset_SN_simulate

end module coop_forecast_mod
