module coop_pertobj_mod
  use coop_wrapper_background
  implicit none
#include "constants.h"

  type coop_pert_hierarchy
     COOP_SHORT_STRING::name
     COOP_INT:: nq = 1
     COOP_INT::m = 0
     COOP_INT::s=0
     COOP_INT::lmax = -1
     COOP_REAL,dimension(:,:),allocatable::var
     COOP_REAL,dimension(:),allocatable::q
   contains
     procedure::init => coop_pert_heirarchy_initialize
     procedure::free => coop_pert_heirarchy_free
  end type coop_pert_hierarchy

contains

  subroutine coop_pert_hierarchy_free(this)
    class(coop_pert_hierarchy)::this
    if(allocated(this%var))deallocate(this%var)
    if(allocated(this%q))deallocate(this%q)
    this%lmax = -1
    this%nq = 1
  end subroutine coop_pert_hierarchy_free

  subroutine coop_pert_hierarchy_initialize(this, lmax, nq, q)
    !!truncate or extend the hierarchy
    COOP_INT, optional::nq
    COOP_REAL,dimension(:), optional::q
    class(coop_pert_hierarchy)::this
    COOP_INT lmax
    COOP_REAL,dimension(:,:),allocatable::tmp
    if(present(nq) )then
       if(.not. present(q)) stop "pert_hierarchy_init: need q for nq parameter"
       call this%free()
       this%nq = nq
       this%lmax = lmax
       allocate(this%var(0:lmax, nq), this%q(nq))
       this%q = q
       this%var = 0.d0
       return
    endif
    if(allocated(this%var))then
       if(this%lmax .eq. lmax)return
       allocate(tmp(0:this%lmax, this%nq))
       tmp = this%var
       deallocate(this%var)
       allocate(this%var(0:lmax, this%nq))
       this%var = 0.d0
       this%var(0: min(lmax, this%lmax),:) = tmp(0: min(lmax, this%lmax),:) 
       this%lmax = lmax
       deallocate(tmp)
    else
       this%lmax = lmax
       allocate(this%var(0:lmax, this%nq))
       this%var = 0.d0
    endif
  end subroutine coop_pert_hierarchy_initialize

  
end module coop_pertobj_mod
