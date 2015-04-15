module coop_forecast_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
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
   contains
     procedure::LogLike =>coop_dataset_SN_loglike
!     procedure::import => coop_dataset_SN_import
!     procedure::export => coop_dataset_SN_export   
  end type coop_DataSet_SN

  
  type coop_MCMC_params
     COOP_INT::proc_id = 0
     type(coop_cosmology_firstorder),pointer::cosmology => null()   
     type(coop_file)::chainfile
     type(coop_file)::logfile
     COOP_STRING::prefix
     COOP_STRING::form
     COOP_INT:: n = 0
     COOP_INT:: fulln = 0
     COOP_REAL:: proposal_length = 1.d0
     COOP_REAL::loglike = coop_LogZero
     COOP_REAL::loglike_proposed = coop_LogZero
     COOP_REAL::mult = 0.d0
     COOP_INT,dimension(:),allocatable::used
     COOP_REAL,dimension(:),allocatable::fullparams
     COOP_REAL,dimension(:),allocatable::params
     COOP_SHORT_STRING, dimension(:), allocatable::name
     COOP_SHORT_STRING, dimension(:), allocatable::tex     
     COOP_REAL,dimension(:),allocatable::lower
     COOP_REAL,dimension(:),allocatable::upper
     COOP_REAL,dimension(:),allocatable::center
     COOP_REAL,dimension(:),allocatable::width     
     COOP_REAL,dimension(:,:),allocatable::covmat
     COOP_REAL,dimension(:,:),allocatable::propose
     COOP_REAL,dimension(:),allocatable::params_saved
     COOP_INT::accept, reject
   contains
     procedure::init => coop_MCMC_params_init
     procedure::MCMC_step => coop_MCMC_params_MCMC_step
  end type coop_MCMC_params



contains


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
       do i=1, this%fulln
          read(fp%unit, *) val
          
       enddo
       call fp%close()
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
    if(allocated(this%fullparams))deallocate(this%used, this%fullparams, this%params, this%lower, this%upper, this%center, this%width, this%covmat, this%propose, this%params_saved, this%name, this%tex)
    allocate(this%fullparams(this%fulln), this%name(this%fulln), this%tex(this%fulln), lower(this%fulln), upper(this%fulln), width(this%fulln), iniwidth(this%fulln), center(this%fulln))
    do i= 1, this%fulln
       this%name(i) = trim(pn%key(i))
       this%tex(i) = trim(pn%val(i))
       call coop_dictionary_lookup(dict, "param["//trim(this%name(i))//"]", val)
       if(trim(val).eq."")then
          write(*,*) "key  param["//trim(this%name(i))//"] is not found in the ini file"//trim(ini)
          stop
       elseif(index(trim(adjustl(val)), " "//coop_tab) .eq. 0)then
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
       else
          this%fullparams(i) = center(i)
       endif
    enddo
    this%n = count(width .ne. 0.d0)    
    allocate(this%used(this%n), this%params(this%n), this%lower(this%n), this%upper(this%n), this%center(this%n), this%width(this%n), this%covmat(this%n, this%n), this%propose(this%n, this%n), this%params_saved(this%n))
    iused= 0
    do i = 1, this%fulln
       if(width(i) .ne. 0.d0)then
          iused = iused + 1
          this%used(iused) = i
       endif
    enddo
    this%center = center(this%used)
    this%width = width(this%used)
    this%upper = upper(this%used)
    this%lower = lower(this%used)

    call coop_dictionary_lookup(dict, "propose_matrix", val)
    if(trim(val).ne."")then
       call fp%open("val", "r")
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
          write(fp%unit, "(A32, 2E16.7)") trim(this%name(i)), lower(i), upper(i)
       enddo
       call fp%close()
    endif
    deallocate(center, lower, upper, width, iniwidth)    
    call dict%free()
    call pn%free()

  end subroutine coop_MCMC_params_init
  

  subroutine coop_MCMC_params_MCMC_step(this, loglike)
    class(coop_MCMC_params)::this
    COOP_REAL loglike, vec(this%n), length
    external loglike
    if(this%accept+this%reject .eq.  0 .and. this%logfile%unit .ne. 0 )then  !!write the settings
       write(this%logfile%unit, "(2I8)") this%n, this%fulln
       write(this%logfile%unit, "("//COOP_STR_OF(this%n)//"I8)") this%used
       write(this%logfile%unit, "("//COOP_STR_OF(this%fulln)//"E16.7)") this%fullparams
    endif
    this%params_saved = this%params
    vec =  coop_random_vector(this%n)*(this%proposal_length*coop_random_Gaussian())       
    this%params = this%params + matmul(this%propose, vec)
    this%fullparams(this%used) = this%params           
    this%loglike_proposed = loglike(this)
    if(this%loglike_proposed - this%loglike .lt. coop_random_exp())then
       this%accept = this%accept + 1
       if(this%chainfile%unit .ne. 0)then
          write(this%chainfile%unit, trim(this%form)) this%mult, this%loglike, this%params_saved
       endif
       this%loglike = this%loglike_proposed
       this%mult  = 1.d0
       return
    else
       this%params = this%params_saved
       this%fullparams(this%used) = this%params       
       this%reject = this%reject + 1
       this%mult = this%mult + 1.d0
    endif
  end subroutine coop_MCMC_params_MCMC_step

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
    COOP_REAL::loglike
    COOP_REAL::mu_theory(this%n), r_theory(this%n), Mbar
    COOP_INT i
    if(this%n .eq. 0 .or. this%off)then
       logLike = 0.d0
       return
    endif
    if(associated(mcmc%cosmology))then
       !$omp parallel do
       do i=1, this%n
          mu_theory(i) = 5.d0*log10(mcmc%cosmology%luminosity_distance(1.d0/(1.d0+this%z(i)))/mcmc%cosmology%H0Mpc())
       enddo
       !$omp end parallel do
    else
       r_theory(1) = coop_integrate(drz, 0.d0, this%z(1))
       !!then I will interpret mcmc%params as (/ Omega_m, Omega_k, w0, wa /)
       do i=2, this%n
          r_theory(i) = r_theory(i-1) + coop_integrate(drz, this%z(i-1), this%z(i))
       enddo
       mu_theory = 5.d0*log10(2997.92458d0*r_theory)       
    endif
    Mbar = sum((mu_theory - this%mu)*this%invdmusq)/this%suminvdmusq
    loglike = sum((mu_theory-this%mu-Mbar)*(mu_theory-this%mu+Mbar)*this%invdmusq)/2.d0
  contains

    function drz(z)
      COOP_REAL z, drz
      drz = 1.d0/sqrt( mcmc%params(1)*(1.d0+z)**3 + mcmc%params(2)*(1.d0+z)**2 + (1.d0-mcmc%params(1))*(1.d0+z)**(3.d0*(1.d0+mcmc%params(3)+mcmc%params(4)))*exp(-3.d0*mcmc%params(4)*z/(1.d0+z)) )
    end function drz
    
  end function coop_dataset_SN_loglike
  
end module coop_forecast_mod
