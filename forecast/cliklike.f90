module coop_clik_mod
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
     logical::initialized = .false.
     COOP_STRING::filename=""
     type(clik_object)::clikid
     logical::is_lensing = .false.
     integer, dimension(6):: has_cl = 0
     integer, dimension(6):: lmax = -1
     integer, dimension(7):: lensing_lmaxs = -1     
     integer::numnames = 0
     integer::n_tot = 0
     character(LEN=256),dimension(:),pointer::names=>null()
     COOP_REAL,dimension(:),allocatable::cl_and_pars, pars
   contains
     procedure::init => coop_clik_object_init
     procedure::free => coop_clik_object_free
     procedure::print_names => coop_clik_object_print_names
     procedure::loglike => coop_clik_object_loglike
     procedure::set_cl_and_pars => coop_clik_object_set_cl_and_pars
  end type coop_clik_object


contains


  subroutine coop_clik_object_init(this, filename)
    class(coop_clik_object)::this
    COOP_UNKNOWN_STRING::filename
    COOP_INT i
    call this%free()
    this%initialized = .true.
    this%filename = trim(filename)    
#ifdef HAS_CLIK
    call clik_try_lensing(this%is_lensing, filename)
    if(this%is_lensing)then
       call clik_lensing_init(this%clikid, filename)
       !!clik setup
       call clik_lensing_get_lmax(this%clikid, this%lmax(1))
       if(this%lmax(1) .eq. -1)then  !!try plik
          call clik_lensing_get_lmaxs(this%clikid, this%lensing_lmaxs)
       else
          this%lensing_lmaxs(1:2) = this%lmax(1)
          this%lensing_lmaxs(3:7) = -1
       endif
       this%numnames = clik_lensing_get_extra_parameter_names(this%clikid, this%names)
       this%n_tot = sum(this%lensing_lmaxs(1:7)+1) + this%numnames
    else
       call clik_init(this%clikid, filename)
       call clik_get_has_cl(this%clikid, this%has_cl)
       call clik_get_lmax(this%clikid, this%lmax)
       this%numnames = clik_get_extra_parameter_names(this%clikid, this%names)
       this%n_tot = sum(this%lmax+1)+this%numnames
    endif
    do i=1, this%numnames
       call coop_convert_to_Fortran_string(this%names(i))
       this%names(i) = trim(adjustl(this%names(i)))
    enddo
    allocate(this%cl_and_pars(this%n_tot), this%pars(this%numnames))       
#else
    stop "CLIK is not installed"
#endif


  end subroutine coop_clik_object_init

  subroutine coop_clik_object_free(this)
    class(coop_clik_object)::this
    if(trim(this%filename).eq."")return
    this%filename=""
    if(associated(this%names))deallocate(this%names)
    this%numnames = 0
    if(allocated(this%cl_and_pars))deallocate(this%cl_and_pars)
    if(allocated(this%pars))deallocate(this%pars)    
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
    this%initialized = .false.
  end subroutine coop_clik_object_free

  !!for lensing Cls(:, 1) = Cl(TT), Cls(:, 2) = Cl(PP)
  !!for 
  subroutine coop_clik_object_set_cl_and_pars(this, Cls, pars)
    class(coop_clik_object)::this
    COOP_INT, parameter::    lmin = 2
    COOP_REAL, dimension(:,lmin:)::Cls
    COOP_REAL, dimension(:), optional::pars
    COOP_INT::istart, i, lmax, ncls, l
    if(.not. this%initialized)return
    if(trim(this%filename) .eq. "") stop "set_cl_and_pars: object not initialized"
    istart = 1
    this%cl_and_pars = 0.d0
    ncls=size(Cls, 1)
    if(ncls .ne. coop_num_cls) stop "set_cl_and_pars: input Cls is different from standard coop Cls format"
    lmax = size(Cls,2)+(lmin-1)
    if(this%is_lensing)then
       if(ncls .lt. 2)then
          stop "set_cl_and_pars: For lensing Cls(1, :) = Cl(TT), Cls(2,:) = Cl(PP)"
       endif

       if(lmax .lt. maxval(this%lensing_lmaxs))then
          write(*,*) trim(this%filename)//" lmax = ",maxval(this%lensing_lmaxs)
          stop "set_cl_and_pars: need to increase lmax"
       endif
       if(lmin .le. this%lensing_lmaxs(1))then
          this%cl_and_pars(istart+lmin:istart+this%lensing_lmaxs(1))= Cls(coop_index_ClLenLen,lmin:this%lensing_lmaxs(1))/(7.430422d12)
       endif
       istart = istart+this%lensing_lmaxs(1)+1                 
       do i = 2, 5
          if(this%lensing_lmaxs(i) .ge. lmin)then
             this%cl_and_pars(istart+lmin:istart+this%lensing_lmaxs(1)) = Cls(i-1, lmin:this%lensing_lmaxs(i))
          endif
          istart = istart+this%lensing_lmaxs(i)+1          
       enddo
    else
       if(lmax .lt. maxval(this%lmax))then
          write(*,*) trim(this%filename)//" lmax = ", maxval(this%lmax)
          stop "set_cl_and_pars: need to increase lmax"
       endif
       do i=1, 4
          if(this%lmax(i) .ge. lmin)then
             this%cl_and_pars(istart + lmin : istart + this%lmax(i)) = Cls(i, lmin : this%lmax(i))
          endif
          istart = istart + this%lmax(i) + 1
       enddo
       istart = istart + this%lmax(5) + this%lmax(6) + 2
    endif
    if(present(pars))then
       if(size(pars).lt. this%numnames)then
          call this%print_names()
          stop "set_cl_and_pars: not enough nuisance parameters"
       endif
       this%cl_and_pars(istart:istart + this%numnames - 1) = pars(1:this%numnames)           
    else
       if(this%numnames .gt. 0)then
          call this%print_names()
          stop "set_cl_and_pars: need to pass the nuisance parameters"
       endif
    endif
  end subroutine coop_clik_object_set_cl_and_pars

  function coop_clik_object_loglike(this) result(loglike)
    class(coop_clik_object)::this
    COOP_REAL loglike
    if(.not. this%initialized)then
       loglike = 0.d0
       return
    endif
#ifdef HAS_CLIK
    if(this%is_lensing)then
       loglike = -clik_lensing_compute(this%clikid, this%cl_and_pars)
    else
       loglike = -clik_compute(this%clikid, this%cl_and_pars)
    endif
    !write(*,*) "check: ", size(this%cl_and_pars), sum(this%cl_and_pars), loglike
#else
    loglike = coop_LogZero
#endif    
  end function coop_clik_object_loglike

  subroutine coop_clik_object_print_names(this, with_param)
    class(coop_clik_object)::this
    COOP_INT ::i
    logical wp
    logical,optional::with_param
    if(.not. this%initialized)return
    if(this%is_lensing)then
       write(*,*) "This is a lensing likelihood"
       write(*,*) "PhiPhi from l=0 to", this%lensing_lmaxs(1)       
       write(*,*) "TT from l=0 to", this%lensing_lmaxs(2)
       write(*,*) "EE from l=0 to", this%lensing_lmaxs(3)
       write(*,*) "BB from l=0 to", this%lensing_lmaxs(4)
       write(*,*) "TE from l=0 to", this%lensing_lmaxs(5)
       write(*,*) "EB from l=0 to", this%lensing_lmaxs(6)
       write(*,*) "TB from l=0 to", this%lensing_lmaxs(7)
    else
       write(*,*) "TT from l=0 to", this%lmax(1)
       write(*,*) "EE from l=0 to", this%lmax(2)
       write(*,*) "BB from l=0 to", this%lmax(3)
       write(*,*) "TE from l=0 to", this%lmax(4)
       write(*,*) "EB from l=0 to", this%lmax(5)
       write(*,*) "TB from l=0 to", this%lmax(6)
    endif
    write(*,*) "=== "//COOP_STR_OF(this%numnames)//" parameters for "//trim(this%filename)//"==="
    if(present(with_param))then
       wp = with_param
    else
       wp = .false.
    endif
    do i=1, this%numnames
       if(wp)then
          write(*,"(A)") "param["//trim(this%names(i))//"] = "
       else
          write(*,"(A)") trim(this%names(i))
       endif
    enddo
    write(*,*) "=== End of parameter list ==="    
  end subroutine coop_clik_object_print_names
  
  
end module coop_clik_mod
