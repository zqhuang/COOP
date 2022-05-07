module coop_gas_mod
  use coop_wrapper_utils
#include "constants.h"

#define SPACE_DIM 2
  COOP_INT,parameter::coop_gas_box_nt = 16
  COOP_INT,parameter::coop_gas_dim = SPACE_DIM

  type  coop_gas_indices
     COOP_INT::n
     COOP_INT::nmax
     COOP_INT,dimension(:),allocatable::i
   contains
     procedure::free => coop_gas_indices_free
     procedure::init => coop_gas_indices_init
     procedure::insert => coop_gas_indices_insert
  end type coop_gas_indices

  type coop_gas_box
     logical::do_collision = .true.
     COOP_INT::n
     COOP_REAL::L, r, twor, fourr2, omega, omega2, Lmr, eps, tsize
     COOP_REAL,dimension(:,:),allocatable::x,v,f
#if SPACE_DIM == 2
     type(coop_gas_indices),dimension(coop_gas_box_nt, coop_gas_box_nt)::ix
#elif SPACE_DIM == 3
     type(coop_gas_indices),dimension(coop_gas_box_nt, coop_gas_box_nt, coop_gas_box_nt)::ix
#endif
   contains
     procedure::free => coop_gas_box_free
     procedure::init => coop_gas_box_init
     procedure::slice => coop_gas_box_slice
     procedure::evolve_x => coop_gas_box_evolve_x
     procedure::evolve_v => coop_gas_box_evolve_v
  end type coop_gas_box

contains 

  subroutine coop_gas_indices_init(this, nmax)
    class(coop_gas_indices)::this
    COOP_INT::nmax
    call this%free()
    allocate(this%i(nmax))
    this%nmax = nmax
    this%n = 0
  end subroutine coop_gas_indices_init

  subroutine coop_gas_indices_free(this)
    class(coop_gas_indices)::this
    this%n = 0
    this%nmax = 0
    COOP_DEALLOC(this%i)
  end subroutine coop_gas_indices_free

  subroutine coop_gas_indices_insert(this, i)
    class(coop_gas_indices)::this
    COOP_INT::i
    COOP_INT,dimension(:),allocatable::tmp
    if(this%n .ge. this%nmax)then
       if(this%nmax .gt. 0)then
          allocate(tmp(this%nmax))
          tmp = this%i
          COOP_DEALLOC(this%i)
          this%nmax = this%nmax + 32
          allocate(this%i(this%nmax))
          this%i(1:this%n) = tmp
          deallocate(tmp)
       else
          this%nmax = 32
          allocate(this%i(this%nmax))
       endif
    endif
    this%n = this%n + 1
    this%i(this%n) = i
  end subroutine coop_gas_indices_insert

  subroutine coop_gas_box_evolve_x(this, dt)
    class(coop_gas_box)::this
    COOP_REAL::dt
    COOP_INT::i, i1, i2, i3
    this%x = this%x + this%v * dt
  end subroutine coop_gas_box_evolve_x


  subroutine coop_gas_box_slice(this)
    class(coop_gas_box)::this
    COOP_INT::i, i1, i2, i3
    this%ix%n = 0
    do i=1, this%n
       i1 = min(max(nint(this%x(1, i)/this%tsize+0.5d0), 1), coop_gas_box_nt)
       i2 = min(max(nint(this%x(2, i)/this%tsize+0.5d0), 1), coop_gas_box_nt)
#if SPACE_DIM == 2
       call this%ix(i1, i2)%insert(i)
#elif SPACE_DIM == 3
       i3 =  min(max(nint(this%x(3, i)/this%tsize+0.5d0), 1), coop_gas_box_nt)
       call this%ix(i1, i2, i3)%insert(i)
#endif       
    enddo

  end subroutine coop_gas_box_slice

  subroutine coop_gas_box_evolve_v(this, dt, pressure)
    class(coop_gas_box)::this
    COOP_REAL::dt
    COOP_REAL::pressure
    COOP_INT::i, dim, i1, i2, i3, j1, j2, j3, k
    COOP_REAL::r2, force(3), r, dp
    pressure = 0.d0
    do dim = 1, SPACE_DIM
       dp = 0.d0
       !$omp parallel do reduction(+:dp)
       do i=1, this%n
          if(this%x(dim, i).gt.this%Lmr)then
             this%f(dim, i) = -this%omega2*(this%x(dim, i) - this%Lmr)
             dp = dp + abs(this%f(dim, i))
          elseif(this%x(dim, i).lt.this%r)then
             this%f(dim, i) = -this%omega2*(this%x(dim, i)-this%r)
             dp = dp + abs(this%f(dim, i))
          else
             this%f(dim, i) = 0.d0
          endif
       enddo
       !$omp end parallel do
       pressure  = pressure + dp
    enddo
    pressure = pressure/(2.d0*SPACE_DIM)/this%L**(SPACE_DIM-1)

    if(this%do_collision)then
       call this%slice()
       do i=1, this%n
          i1 = min(max(nint(this%x(1, i)/this%tsize+0.5d0), 1), coop_gas_box_nt)
          i2 = min(max(nint(this%x(2, i)/this%tsize+0.5d0), 1), coop_gas_box_nt)
#if SPACE_DIM == 2
          do j2 = max(1, i2-1), min(coop_gas_box_nt, i2+1)
             do j1 = max(1, i1-1), min(coop_gas_box_nt, i1+1)
                do k = 1, this%ix(j1, j2)%n
                   if(this%ix(j1, j2)%i(k).ne. i)then
                      r2 = sum((this%x(:, i) - this%x(:, this%ix(j1, j2)%i(k)))**2)
                      if(r2 .lt. this%fourr2)then
                         r = sqrt(r2)
                         this%f(:, i) = this%f(:, i) + this%omega2 * (this%twor - r) * (this%x(:, i) - this%x(:, this%ix(j1, j2)%i(k)))/(r+this%eps)
                      endif
                   endif
                enddo
             enddo
          enddo
#elif SPACE_DIM == 3
          i3  = min(max(nint(this%x(3, i)/this%tsize+0.5d0), 1), coop_gas_box_nt)
          do j3 = max(1, i3-1), min(coop_gas_box_nt, i3+1)
             do j2 = max(1, i2-1), min(coop_gas_box_nt, i2+1)
                do j1 = max(1, i1-1), min(coop_gas_box_nt, i1+1)
                   do k = 1, this%ix(j1, j2, j3)%n
                      if(this%ix(j1, j2, j3)%i(k).ne. i)then
                         r2 = sum((this%x(:, i) - this%x(:, this%ix(j1, j2, j3)%i(k)))**2)
                         if(r2 .lt. this%fourr2)then
                            r = sqrt(r2)
                            this%f(:, i) = this%f(:, i) + this%omega2 * (this%twor - r) * (this%x(:, i) - this%x(:, this%ix(j1, j2, j3)%i(k)))/(r+this%eps)
                         endif
                      endif
                   enddo
                enddo
             enddo
          enddo
#endif
       enddo
    endif
    this%v = this%v + this%f*dt
  end subroutine coop_gas_box_evolve_v

  
  subroutine coop_gas_box_free(this)
    class(coop_gas_box)::this
    COOP_DEALLOC(this%x)
    COOP_DEALLOC(this%v)
    COOP_DEALLOC(this%f)
    this%n  = 0
  end subroutine coop_gas_box_free

  subroutine coop_gas_box_init(this, n, L, r)
    class(coop_gas_box)::this
    COOP_REAL::L, r
    COOP_INT::n, i1, i2, i3, nmax
    call this%free()
    this%n   = n
    this%L = L
    this%r = r
    this%twor = 2.d0*r
    this%fourr2 = this%twor**2
    this%Lmr = this%L - this%r
    this%eps = this%r/1.d2
    this%omega = 10.d0/this%r
    this%omega2 = this%omega**2
    this%tsize = this%L/coop_gas_box_nt
    allocate(this%x(SPACE_DIM, n), this%v(SPACE_DIM, n), this%f(SPACE_DIM, n))
    if(this%do_collision)then
       nmax = (this%n / (coop_gas_box_nt**SPACE_DIM) + 1)*4
       do i1=1, coop_gas_box_nt
          do i2=1, coop_gas_box_nt
#if SPACE_DIM == 2
             call this%ix(i2, i1)%init(nmax)
#elif SPACE_DIM == 3
             do i3= 1, coop_gas_box_nt
                call this%ix(i3, i2, i1)%init(nmax)
             enddo
#endif
          enddo
       enddo
    endif
  end subroutine coop_gas_box_init
  

#undef SPACE_DIM
end module coop_gas_mod
