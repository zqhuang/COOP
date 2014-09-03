module coop_pertobj_mod
  use coop_wrapper_background
  implicit none
#include "constants.h"

  COOP_INT, parameter:: coop_pert_default_lmax = 32  
  COOP_INT, parameter:: coop_pert_default_mmax = 2
  COOP_INT, parameter:: coop_pert_default_smax = 2
  COOP_INT, parameter:: coop_pert_default_nq = 5
  COOP_REAL, dimension(coop_pert_default_nq),parameter:: coop_pert_default_q = coop_fermion_int_q5
  COOP_REAL, dimension(coop_pert_default_nq),parameter:: coop_pert_default_q_kernel = coop_fermion_int_kernel5 


  COOP_INT, parameter::coop_pert_imetric = 1
  COOP_INT, parameter::coop_pert_ibaryon = 2
  COOP_INT, parameter::coop_pert_icdm = 3
  COOP_INT, parameter::coop_pert_iT = 4
  COOP_INT, parameter::coop_pert_iE = 5
  COOP_INT, parameter::coop_pert_iB = 6
  COOP_INT, parameter::coop_pert_iNu = 7
  COOP_INT, parameter::coop_pert_imassiveNu = 8
  COOP_INT, parameter::coop_pert_ide = 9
  COOP_INT, parameter::coop_pert_nspecies = coop_pert_ide


  type coop_pert_hierarchy
     COOP_INT:: nq = 1
     COOP_INT::m = 0
     COOP_INT::s=0
     COOP_INT::lmax = -1
     COOP_INT::lmin_used = -1
     COOP_INT::lmax_used = -1
     COOP_INT::nvars_used = 0
     COOP_INT::nvars_used_per_q = 0
     COOP_REAL,dimension(:,:),allocatable::var
     COOP_REAL,dimension(:),allocatable::q
   contains
     procedure::init => coop_pert_hierarchy_initialize
     procedure::nvars => coop_pert_hierarchy_nvars
     procedure::free => coop_pert_hierarchy_free
  end type coop_pert_hierarchy


  type coop_standard_o1pert
     COOP_SHORT_STRING::initial_conditions = "adiabatic"
     type(coop_pert_hierarchy),dimension(coop_pert_nspecies)::species
     COOP_INT, dimension(coop_pert_nspecies)::istart, iend
     COOP_REAL::k, phi, pi, phidot, pidot, vdot
     logical:: tight_coupling = .true.
     logical:: has_massiveNu = .false.
     COOP_INT :: de_pert_model = 0
     logical:: late_approx = .false.
     COOP_INT :: m = 0
     COOP_REAL,dimension(:),allocatable::y
     COOP_REAL,dimension(:, :),allocatable::ode_w
     COOP_REAL::ode_c(24)
     COOP_INT::ode_ind, ode_nvars
   contains
     procedure::init => coop_standard_o1pert_initialize
     procedure::set_ode_index => coop_standard_o1pert_set_ode_index
     procedure::read_var => coop_standard_o1pert_read_var
     procedure::write_var => coop_standard_o1pert_write_var
     procedure::free => coop_standard_o1pert_free
  end type coop_standard_o1pert

contains

  subroutine coop_standard_o1pert_initialize(this, m)
    class(coop_standard_o1pert)::this
    COOP_INT::m, i
    logical::tight_coupling
    this%m = m
    this%species%m = m
    this%species(coop_pert_iE)%s = 2  !!override the default spin 0
    this%species(coop_pert_iB)%s = 2
    do i=1, coop_pert_nspecies
       select case(i)
       case(coop_pert_imetric)
          call this%species(i)%init(lmax = 4)          
       case(coop_pert_ibaryon, coop_pert_icdm)
          call this%species(i)%init(lmax = 2)
       case(coop_pert_imassiveNu)
          call this%species(i)%init(lmax = coop_pert_default_lmax, nq =  coop_pert_default_nq, q = coop_pert_default_q )
       case(coop_pert_ide)
          call this%species(i)%init(lmax = 4)
       case default
          call this%species(i)%init(lmax = coop_pert_default_lmax)
       end select
    enddo
  end subroutine coop_standard_o1pert_initialize


  subroutine coop_standard_o1pert_free(this)
    class(coop_standard_o1pert)::this
    COOP_INT i
    do i = 1, coop_pert_nspecies
       call this%species(i)%free()
    enddo
  end subroutine coop_standard_o1pert_free


  subroutine coop_standard_o1pert_set_ode_index(this)
    class(coop_standard_o1pert)::this
    COOP_INT i
    select case(this%m)
    case(0,1)
       this%species(coop_pert_imetric)%lmax_used = 1
    case(2)
       this%species(coop_pert_imetric)%lmax_used = 3
    case default
       call coop_return_error("set_ode_index", "invalid m: "//trim(coop_num2str(this%m)), "stop")
    end select
       
    this%species(coop_pert_icdm)%lmax_used = 1
    this%species(coop_pert_iNu)%lmax_used = 12
    this%species(coop_pert_ibaryon)%lmax_used = 1
    if(this%tight_coupling)then
       this%species(coop_pert_iT)%lmax_used = 0
       this%species(coop_pert_iE)%lmax_used = 0
       this%species(coop_pert_iB)%lmax_used = 0
    else
       if(this%late_approx)then
          this%species(coop_pert_iT)%lmax_used = -1
          this%species(coop_pert_iE)%lmax_used = -1
          this%species(coop_pert_iB)%lmax_used = -1
       else
          this%species(coop_pert_iT)%lmax_used = 16
          this%species(coop_pert_iE)%lmax_used = 12
          if(this%m .gt. 0)then
             this%species(coop_pert_iB)%lmax_used =  this%species(coop_pert_iE)%lmax_used
          else
             this%species(coop_pert_iB)%lmax_used = -1
          endif
       endif
    endif
    if(this%has_massiveNu)then
       this%species(coop_pert_imassiveNu)%lmax_used = 10
    else
       this%species(coop_pert_imassiveNu)%lmax_used = -1
    endif
    select case(this%de_pert_model)
    case(COOP_DE_PERT_NONE)
       this%species(coop_pert_ide)%lmax_used = -1
    case(COOP_DE_PERT_FLUID)
       this%species(coop_pert_ide)%lmax_used = 1
    case(COOP_DE_PERT_PPF)
       this%species(coop_pert_ide)%lmax_used = 0
    case default
       call coop_return_error("set_ode_index", "unknown de_pet_model: "//trim(coop_num2str(this%de_pert_model)), "stop")
    end select
    this%istart(1) = 1
    this%iend(1) = this%istart(1) + this%species(1)%nvars() - 1
    do i=2, coop_pert_nspecies
       this%istart(i) = this%iend(i-1)+1
       this%iend(i) = this%istart(i) - 1 + this%species(i)%nvars()
    enddo
    this%ode_nvars = this%iend(coop_pert_nspecies)
    if(allocated(this%ode_w))deallocate(this%ode_w)
    if(allocated(this%y))deallocate(this%y)
    allocate(this%ode_w(this%ode_nvars, 9), this%y(this%ode_nvars))
    this%ode_ind = 1
    this%ode_c = 0.
    this%ode_w = 0
  end subroutine coop_standard_o1pert_set_ode_index

  subroutine coop_standard_o1pert_read_var(this, var)
    COOP_REAL,dimension(:),optional::var
    class(coop_standard_o1pert)::this
    COOP_INT i
    if(present(var))then
       do i=1, coop_pert_nspecies
          if(this%species(i)%nvars_used .gt. 0)then
             this%species(i)%var(this%species(i)%lmin_used:this%species(i)%lmax_used, :) = reshape(var(this%istart(i):this%iend(i)), (/ this%species(i)%nvars_used_per_q, this%species(i)%nq /))
          endif
       enddo
    else
       do i=1, coop_pert_nspecies
          if(this%species(i)%nvars_used .gt. 0)then
             this%species(i)%var(this%species(i)%lmin_used:this%species(i)%lmax_used, :) = reshape(this%y(this%istart(i):this%iend(i)), (/ this%species(i)%nvars_used_per_q, this%species(i)%nq /))
          endif
       enddo
    endif
  end subroutine coop_standard_o1pert_read_var

  subroutine coop_standard_o1pert_write_var(this, var)
    class(coop_standard_o1pert)::this
    COOP_REAL,dimension(:), optional::var
    COOP_INT i
    if(present(var))then
       do i=1, coop_pert_nspecies
          if(this%species(i)%nvars_used .gt. 0)then
             var(this%istart(i):this%iend(i)) = reshape(this%species(i)%var(this%species(i)%lmin_used:this%species(i)%lmax_used, :), (/ this%species(i)%nvars_used /) )
          endif
       enddo
    else
       do i=1, coop_pert_nspecies
          if(this%species(i)%nvars_used .gt. 0)then
             this%y(this%istart(i):this%iend(i)) = reshape(this%species(i)%var(this%species(i)%lmin_used:this%species(i)%lmax_used, :), (/ this%species(i)%nvars_used /) )
          endif
       enddo
    endif
  end subroutine coop_standard_o1pert_write_var

  subroutine coop_pert_hierarchy_free(this)
    class(coop_pert_hierarchy)::this
    if(allocated(this%var))deallocate(this%var)
    if(allocated(this%q))deallocate(this%q)
    this%lmax = -1
    this%lmax_used = -1
    this%lmin_used = -1
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
       if(lmax .ge. 0)   allocate(this%var(-1:lmax, nq))
       if(nq.ge.1)then
          allocate(this%q(nq))
          this%q = q
          this%var = 0.d0
       endif
       this%lmax_used = this%lmax - 1
       this%lmin_used = max(this%m, this%s)
       return
    else
       this%nq = 1
       if(.not. allocated(this%q)) allocate(this%q(1))
       this%q = 1.d0
    endif
    if(allocated(this%var))then
       if(this%lmax .eq. lmax)return
       allocate(tmp(-1:this%lmax, this%nq))
       tmp = this%var
       deallocate(this%var)
       if(lmax .ge. 0)then
          allocate(this%var(-1:lmax, this%nq))
          this%var = 0.d0
          this%var(-1: min(lmax, this%lmax),:) = tmp(-1: min(lmax, this%lmax),:) 
       endif
       this%lmax = lmax
       deallocate(tmp)
    else
       this%lmax = lmax
       if(lmax .ge. 0)then
          allocate(this%var(-1:lmax, this%nq))
          this%var = 0.d0
       endif
    endif
    this%lmax_used = this%lmax - 1
    this%lmin_used = max(this%m, this%s)
  end subroutine coop_pert_hierarchy_initialize

  function coop_pert_hierarchy_nvars(this) result(n)
    COOP_INT n
    class(coop_pert_hierarchy)::this
    this%lmin_used = max(this%m, this%s, this%lmin_used)
    this%lmax_used = max(min(this%lmax_used, this%lmax-1), this%lmin_used - 1)
    this%nvars_used_per_q = (this%lmax_used- this%lmin_used+1)
    this%nvars_used = this%nvars_used_per_q * this%nq
    n = this%nvars_used 
  end function coop_pert_hierarchy_nvars


end module coop_pertobj_mod
