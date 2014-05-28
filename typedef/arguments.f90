module coop_type_arguments
  use coop_constants
  use coop_basicutils
  implicit none
#include "constants.h"

private

public::coop_arguments

  type coop_arguments
     COOP_INT n_int, n_real, n_logical
     COOP_INT,dimension(:),allocatable::i
     COOP_REAL,dimension(:),allocatable::r
     logical,dimension(:),allocatable::l
   contains
     procedure::init => coop_arguments_initialize
     procedure::free => coop_arguments_free
  end type coop_arguments

  interface coop_arguments
     procedure coop_arguments_constructor
  end interface coop_arguments

contains


  function coop_arguments_constructor(i, r , l) result(this)
    type(coop_arguments) this
    COOP_INT, dimension(:),optional::i
    COOP_REAL, dimension(:),optional::r
    logical, dimension(:),optional::l
    call this%free
    if(present(i))then
       this%n_int = size(i)
       allocate(this%i(this%n_int))
       this%i = i
    endif
    if(present(r))then
       this%n_real = size(r)
       allocate(this%r(this%n_real))
       this%r = r
    endif
    if(present(l))then
       this%n_logical = size(l)
       allocate(this%l(this%n_logical))
       this%l = l
    endif
  end function coop_arguments_constructor

  subroutine coop_arguments_initialize(this,i,r,l)
    class(coop_arguments) this
    COOP_INT, dimension(:),optional::i
    COOP_REAL, dimension(:),optional::r
    logical, dimension(:),optional::l
    call this%free
    if(present(i))then
       this%n_int = size(i)
       allocate(this%i(this%n_int))
       this%i = i
    endif
    if(present(r))then
       this%n_real = size(r)
       allocate(this%r(this%n_real))
       this%r = r
    endif
    if(present(l))then
       this%n_logical = size(l)
       allocate(this%l(this%n_logical))
       this%l = l
    endif
  end subroutine coop_arguments_initialize
    
  subroutine coop_arguments_free(this)
    class(coop_arguments) this
    if(allocated(this%i))deallocate(this%i)
    if(allocated(this%r))deallocate(this%r)
    if(allocated(this%l))deallocate(this%l)
    this%n_int = 0
    this%n_real = 0
    this%n_logical = 0
  end subroutine coop_arguments_free

end module coop_type_arguments
