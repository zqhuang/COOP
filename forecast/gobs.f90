module coop_gobs_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

 ! COOP_REAL,parameter::hrd_fiducial = 0.7*147.78  !!h * r_drag in Mpc, as measured by Planck
  
  type coop_gobs_object
     COOP_INT::n = 0
     logical::is_diag = .false.
     logical::is_log = .false.
     COOP_STRING,dimension(:),allocatable::name
     COOP_REAL,dimension(:),allocatable::val, z, zs, log_shift
     COOP_REAL,dimension(:),allocatable::val_theory, sigma
     COOP_REAL, dimension(:,:), allocatable::invcov    
   contains
     procedure::loglike => coop_gobs_loglike
     procedure::free => coop_gobs_free
     procedure::init => coop_gobs_init
  end type coop_gobs_object

contains


  subroutine coop_gobs_free(this)
    class(coop_gobs_object)::this
    COOP_DEALLOC(this%name)
    COOP_DEALLOC(this%val)
    COOP_DEALLOC(this%z)
    COOP_DEALLOC(this%log_shift)    
    COOP_DEALLOC(this%zs)
    COOP_DEALLOC(this%sigma)        
    COOP_DEALLOC(this%val_theory)
    COOP_DEALLOC(this%invcov)    
    this%n = 0
    this%is_diag = .false.
    this%is_log = .false.
  end subroutine coop_gobs_free

  subroutine coop_gobs_init(this, filename)
    class(coop_gobs_object)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file)::fp
    COOP_INT::i
    COOP_STRING::line
    type(coop_list_string)::ls
    COOP_STRING::cov_format
    call this%free()
    call fp%open(COOP_DATAPATH(filename))
    write(*,*) "loading data from "//trim(COOP_DATAPATH(filename))
    read(fp%unit, "(A)", END=100, ERR=100) line  !!head information
    write(*,*) trim(line)
    read(fp%unit, *, END=100, ERR=100) this%n
    if(this%n .le. 0) goto 100    
    write(*,*) this%n, " data points"    
    read(fp%unit, *, END=100, ERR=100) this%is_diag
    if(this%is_diag) write(*,*) "covariance is diagonal"
    read(fp%unit, *, END=100, ERR=100) this%is_log
    allocate(this%name(this%n), this%log_shift(this%n), this%val(this%n), this%sigma(this%n), this%val_theory(this%n), this%z(this%n), this%zs(this%n), this%invcov(this%n, this%n))
    this%sigma = 0.d0
    if(this%is_log)then
       read(fp%unit,*, END=100, ERR=100) this%log_shift
    endif
    do i = 1, this%n
       read(fp%unit, "(A)", ERR=100, END=100) line
       call coop_string_to_list(line, ls, deliminator=",", ignore_repeat = .false.)
       this%name(i) = trim(coop_str_NumUpperAlpha(ls%element(1)))
       call ls%get_element_as_real(2, this%val(i))
       if(trim(ls%element(3)) .ne. "") &
            call ls%get_element_as_real(3,this%sigma(i))
       if(trim(ls%element(4)) .ne. "") &
            call ls%get_element_as_real(4, this%z(i))
       if(trim(ls%element(5)) .ne. "") &
            call ls%get_element_as_real(5,this%zs(i))
    enddo
    if(this%is_diag)then  !need to read the covariant matrix
       if(any(this%sigma .le. 0.d0)) goto 100
    else
       read(fp%unit, "(A)", END=100, ERR=100) cov_format
       call coop_str2upper(cov_format)
       do i=1, this%n
          read(fp%unit, *, ERR=100, END=100) this%invcov(i,:)
       enddo
       select case(trim(cov_format))
       case("COV")
          call coop_matrix_inverse(this%invcov)
       case("CORR")
          do i=1, this%n
             this%invcov(i,:) = this%invcov(i,:)*this%sigma(i)
             this%invcov(:,i) = this%invcov(:,i)*this%sigma(i)             
          enddo
          call coop_matrix_inverse(this%invcov)          
       case("INVCOV")
          !do nothing
       case default
          goto 100
       end select
    endif
    call fp%close()
    write(*,*) "data loaded!"
    return
100 write(*,*) "Error in data file: "//trim(filename)
    stop
  end subroutine coop_gobs_init

  function coop_gobs_loglike(this, cosmology) result(loglike)  !!return  - ln(likelihood) = chi^2/2
    COOP_REAL,parameter::eps = 1.e-11
    class(coop_gobs_object)::this
    COOP_REAL::loglike
    type(coop_cosmology_firstorder),optional::cosmology
    COOP_INT::i
    do i=1, this%n
       this%val_theory(i) =  cosmology%general_variable(this%name(i), this%z(i), this%zs(i), cosmology%r_drag)
    enddo
!!$    print*, this%val_theory
!!$    print*, this%val
!!$    print*
!!$    do i=1, size(this%invcov,1)
!!$       print*, this%invcov(i,:)
!!$    enddo
!!$    stop
    if(this%is_log)then
       if(any(this%val_theory .le. this%log_shift + eps))then
          loglike = coop_logzero
          return
       endif
       this%val_theory = log(this%val_theory-this%log_shift)
    endif
    if(this%is_diag)then
       loglike = sum(((this%val - this%val_theory)/this%sigma)**2)/2.0       
    else
       loglike = dot_product(this%val - this%val_theory, matmul(this%invcov, this%val-this%val_theory))/2.0
    endif
    if(this%is_log)then
       loglike = loglike +  sum(this%val_theory)
    endif
    
  end function coop_gobs_loglike
  
end module coop_gobs_mod
