module coop_clik
#ifdef HAS_CLIK  
  use clik
#endif  
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

#ifdef HAS_CLIK
#else  
  !!define the clik objects
  type clik_object
     COOP_STRING::name = "CLIK_NOT_FOUND"
  end type clik_object
#endif  

  type coop_clik_object
     COOP_STRING::filename=""
     type(clik_object)::clikid
     logical::is_lensing = .false.
     integer, dimension(6):: has_cl = 0
     integer, dimension(6):: lmax = -1
     integer, allocatable:: clik_index(:,:)
     integer::numnames = 0
     integer::n_tot = 0
     character(LEN=256),dimension(:),pointer::names
     COOP_REAL,dimension(:),allocatable::cl_and_pars
   contains
     procedure::init => coop_clik_object_init
     procedure::free => coop_clik_object_free
     procedure::loglike => coop_clik_object_loglike
     procedure::set_cl_and_pars => coop_clik_object_set_cl_and_pars
  end type coop_clik_object


contains


  subroutine coop_clik_object_init(this, filename)
    class(coop_clik_object)::this
    COOP_UNKNOWN_STRING::filename
    this%filename = trim(filename)
    call this%free()
#ifdef HAS_CLIK
    call clik_try_lensing(this%is_lensing, filename)
    if(this%is_lensing)then
       call clik_lensing_init(this%clikid, filename)
       call clik_get_lmax(this%clikid, this%lmax(1))
       this%numnames = clik_get_extra_parameter_names(this%clikid, this%names)
       this%n_tot = 2*(this%lmax(1)+1) + this%numnames
    else
       call clik_init(this%clikid, filename)
       call clik_get_has_cl(this%clikid, this%has_cl)
       call clik_get_lmax(this%clikid, this%lmax)
       this%numnames = clik_get_extra_parameter_names(this%clikid, this%names)
       this%n_tot = sum(this%lmax)+6+this%numnames
    endif
#else
    stop "CLIK is not installed"
#endif
    allocate(this%cl_and_pars(this%n_tot))
  end subroutine coop_clik_object_init

  subroutine coop_clik_object_free(this)
    class(coop_clik_object)::this
    this%filename=""
    if(associated(this%names))deallocate(this%names)
    this%numnames = 0
    if(allocated(this%cl_and_pars))deallocate(this%cl_and_pars)
    this%lmax = -1
    this%n_tot = 0
    this%has_cl = 0
#ifdef HAS_CLIK
    if(this%is_lensing)then
       call clik_lensing_cleanup(this%clikid)
    else
       call clik_cleanup(this%clikid)
    endif
#endif    
  end subroutine coop_clik_object_free

  !!for lensing Cls(:, 1) = Cl(TT), Cls(:, 2) = Cl(PP)
  !!for 
  subroutine coop_clik_object_set_cl_and_pars(this, Cls, pars)
    class(coop_clik_object)::this
    COOP_REAL, dimension(:,:)::Cls
    COOP_REAL, dimension(:)::pars
    COOP_INT::istart, i, lmin, lmax, ncls
    if(trim(this%filename) .eq. "") stop "set_cl_and_pars: object not initialized"
    istart = 1
    this%cl_and_pars = 0.d0
    ncls=size(Cls, 2)
    lmin = lbound(Cls,2)
    lmax = ubound(Cls,2)
    if(lmin .gt. 2)stop "set_cl_and_pars: lmin > 2?"    
    if(this%is_lensing)then
       if(ncls .ne. 2)then
          stop "set_cl_and_pars: For lensing Cls(:, 1) = Cl(TT), Cls(:, 2) = Cl(PP)"
       endif

       if(lmax .lt. this%lmax(1))then
          write(*,*) trim(this%filename)//" lmax = ", this%lmax(1)
          stop "set_cl_and_pars: need to increase lmax"
       endif
       if(this%lmax(1) .ge. lmin)then
          this%cl_and_pars(istart+lmin:istart+this%lmax(1)) = Cls(lmin:this%lmax(1), 2)
          istart = istart+this%lmax(1)+1
          this%cl_and_pars(istart+lmin:istart+this%lmax(1)) = Cls(lmin:this%lmax(1), 1)
       endif
    else
       if(ncls .ne. 4 .and. ncls .ne. 6)then
          stop "set_cl_and_pars: For lensing Cls(:, 1) = Cl(TT), Cls(:, 2) = Cl(EE), Cls(:,3) = Cl(BB), Cls(:,4) = Cl(TE) [optional Cls(:,5) = Cl(TB), Cls(:,6) = Cl(EB)]"
       endif
       if(lmin .gt. 2)stop "set_cl_and_pars: lmin is at most 2"
       if(lmax .lt. maxval(this%lmax))then
          write(*,*) trim(this%filename)//" lmax = ", maxval(this%lmax)
          stop "set_cl_and_pars: need to increase lmax"
       endif
       do i = 1, ncls
          if(this%lmax(i) .ge. lmin)then
             this%cl_and_pars(istart + lmin : istart + this%lmax(i)) = Cls(lmin : this%lmax(i), i)
          endif
          istart = istart + this%lmax(i) + 1
          this%cl_and_pars(istart + lmin : istart + this%lmax(1)) = Cls(lmin : this%lmax(1), 1)
       enddo
          
    endif
  end subroutine coop_clik_object_set_cl_and_pars

  function coop_clik_object_loglike(this) result(loglike)
    class(coop_clik_object)::this
    COOP_REAL loglike
#ifdef HAS_CLIK
    if(this%is_lensing)then
       loglike = clik_lensing_compute(this%clikid, this%cl_and_pars)
    else
       loglike = clik_compute(this%clikid, this%cl_and_pars)
    endif
#else
    loglike = coop_LogZero
#endif    
  end function coop_clik_object_loglike
  
end module coop_clik
