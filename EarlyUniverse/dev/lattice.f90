module coop_lattice_mod
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  !!This is HLattice v3.0: a fine-grained lattice for field and coarse-grained lattice for gravity is used.
#include "lattice.h"  

  private

  public:: coop_lattice
  
  type coop_lattice
     COOP_INT::n_H_terms = 3
     COOP_INT::ode_order = 6
     COOP_INT::n = 0
     COOP_INT::nc = 0
     COOP_INT::cs = 0
     COOP_INT::cs3 = 0
     COOP_INT::dim_f = 0
     COOP_REAL,dimension(:,:,:,:),allocatable::f, pi_f
     COOP_REAL,dimension(:,:,:,:),allocatable::beta, pi_beta, hup
     COOP_REAL::dx = 1.d0
     COOP_REAL::y = 1.d0
     COOP_REAL::pi_y = 0.d0
     COOP_REAL::mu = 1.d0  !!d (cosmological time)  = a^mu d (program time)
     COOP_REAL::Mpsq = 1024.d0**2.d0  !!reduced Planck mass^2, in principle you can set it to be any positive number. Practically a proper choice of Mpsq can make the numbers not too big/too small.
   contains
     procedure::free => coop_lattice_free
     procedure::alloc => coop_lattice_alloc
     procedure::wrap_f => coop_lattice_wrap_f
     procedure::wrap_beta => coop_lattice_wrap_beta
     procedure::get_hup => coop_lattice_get_hup
     procedure::wrap_all => coop_lattice_wrap_all
     procedure::evolve => coop_lattice_evolve
     procedure::symp2 => coop_lattice_symp2
     procedure::symp4 => coop_lattice_symp4
     procedure::symp6 => coop_lattice_symp6
     procedure::symp_o2step => coop_lattice_symp_o2step
     procedure::H_split => coop_lattice_H_split
     procedure::H_K1 => coop_lattice_H_K1
     procedure::H_K2 => coop_lattice_H_K2
     procedure::H_K3 => coop_lattice_H_K3     
     procedure::H_P => coop_lattice_H_P
     procedure::K1 => coop_lattice_K1
     procedure::K2 => coop_lattice_K2
     procedure::K3 => coop_lattice_K3
     procedure::P => coop_lattice_P
     procedure::V_f => coop_lattice_V_f
     procedure::G_f => coop_lattice_G_f
     procedure::G_g => coop_lattice_G_g     
  end type coop_lattice
  
contains


  !!============ editable section: define your model here ===========
#define PHI phi(1)
#define CHI phi(2)  
  function coop_lattice_potential(phi) result(V)
    COOP_REAL::phi(:), V
    COOP_REAL, parameter::lambda = 1.d-13 !!lambda
    COOP_REAL, parameter::g2byl2 = 2.d0 !!g^2/lambda^2
    V =  lambda * ( (1.d0/4.d0) * PHI **2 + (g2byl2/2.d0) * CHI **2 ) * PHI**2
  end function coop_lattice_potential

  function coop_lattice_dVdphi(phi) result(dVdphi)
    COOP_REAL::phi(:), dVdphi(size(phi))
    COOP_REAL, parameter::lambda = 1.d-13 !!lambda
    COOP_REAL, parameter::g2byl2 = 2.d0 !!g^2/lambda^2
    dVdphi(1) = lambda*(PHI**2 + g2byl2 * CHI**2) * PHI
    dVdphi(2) = (lambda*g2byl2) * CHI * PHI**2
  end function coop_lattice_dVdphi
  
#undef PHI
#undef CHI  
  !!===========================end of editable section ==============
  

  subroutine coop_lattice_free(this)
    class(coop_lattice)::this    
    if(allocated(this%f))deallocate(this%f)
    if(allocated(this%pi_f))deallocate(this%pi_f)
    if(allocated(this%beta))deallocate(this%beta)
    if(allocated(this%pi_beta))deallocate(this%pi_beta)
    this%n = 0
    this%nc  = 0
    this%dim_f = 0
    this%cs = 0
    this%cs3 = 0
    this%y = 1.d0
    this%pi_y = 0.d0
  end subroutine coop_lattice_free

  subroutine coop_lattice_alloc(this,dim_f, n, nc)
    class(coop_lattice)::this
    COOP_INT::dim_f, n, nc
    if(n .le. nc .or. nc .le. 0 )then
       stop "coop_lattice_alloc failed: field grids must be finer than metric grids"
    elseif(mod(n, nc) .ne. 0)then
       stop "coop_lattice_alloc failed: n/nc must be an integer"
    endif
    call this%free()
    this%dim_f = dim_f
    this%n = n
    this%nc = nc
    this%cs = this%n / this%nc
    this%cs3 = this%cs ** 3
    allocate(this%f(dim_f, -1:n, -1:n, -1:n))
    allocate(this%pi_f(dim_f, 0:n-1, 0:n-1, 0:n-1))
    allocate(this%beta(6, -1:nc, -1:nc, -1:nc))
    allocate(this%hup(6, -1:nc, -1:nc, -1:nc))    
    allocate(this%pi_beta(6, 0:nc-1, 0:nc-1, 0:nc-1))
  end subroutine coop_lattice_alloc

  subroutine coop_lattice_wrap_f(this)
    class(coop_lattice)::this    
    this%f(:, -1, :, :) = this%f(:, this%n - 1, :, :)
    this%f(:, :, -1, :) = this%f(:, :, this%n - 1, :)
    this%f(:, :, :, -1) = this%f(:, :, :, this%n - 1)
    this%f(:, this%n, :, :) = this%f(:, 0, :, :)
    this%f(:, :, this%n, :) = this%f(:, :, 0, :)
    this%f(:, :, :, this%n) = this%f(:, :, :, 0)
  end subroutine coop_lattice_wrap_f


  subroutine coop_lattice_wrap_beta(this)
    class(coop_lattice)::this    
    this%beta(:, -1, :, :) = this%beta(:, this%nc - 1, :, :)
    this%beta(:, :, -1, :) = this%beta(:, :, this%nc - 1, :)
    this%beta(:, :, :, -1) = this%beta(:, :, :, this%nc - 1)
    this%beta(:, this%nc, :, :) = this%beta(:, 0, :, :)
    this%beta(:, :, this%nc, :) = this%beta(:, :, 0, :)
    this%beta(:, :, :, this%nc) = this%beta(:, :, :, 0)
  end subroutine coop_lattice_wrap_beta

  subroutine coop_lattice_wrap_all(this)
    class(coop_lattice)::this    
    call this%wrap_f()
    call this%wrap_beta()
  end subroutine coop_lattice_wrap_all
  
  

  subroutine coop_lattice_evolve(this, dt, nsteps)
    class(coop_lattice)::this
    COOP_REAL::dt
    COOP_INT nsteps
    select case(this%ode_order)
    case(2)
       call this%symp2(dt, nsteps)
    case(4)
       call this%symp4(dt, nsteps)
    case(6)
       call this%symp6(dt, nsteps)
    case default
       stop "coop_lattice: ode_order must be 2, 4, or 6"
    end select
  end subroutine coop_lattice_evolve

  subroutine coop_lattice_symp_o2step(this, dt,c1,c2)
    class(coop_lattice)::this
    COOP_REAL dt,c1,c2
    COOP_INT::i
    do i=2, this%n_H_terms-1
       call this%H_split(c1*dt/2.d0,i)
    enddo
    call this%H_Split(c1*dt, this%n_H_terms)
    do i=this%n_H_terms-1,2,-1
       call this%H_Split(c1*dt/2.d0, i)
    enddo
    call this%H_Split((c1+c2)*dt/2.d0,1)
    return
  end subroutine coop_lattice_symp_o2step

  subroutine coop_lattice_symp2(this, dt,nsteps)
    class(coop_lattice)::this
    COOP_REAL::dt
    COOP_INT::nsteps, j
    call this%H_split(dt/2.d0, 1)
    do j=1,nsteps-1
       call this%symp_o2step(dt, 1.d0, 1.d0)
    enddo
    call this%symp_o2step(dt, 1.d0, 0.d0)
  end subroutine coop_lattice_symp2

  subroutine coop_lattice_symp4(this, dt,nsteps)
    class(coop_lattice)::this
    COOP_REAL::dt
    COOP_INT::nsteps, j
    COOP_REAL,parameter:: c1 = 1.d0/(2.d0 - 2.d0**(1.d0/3.d0))
    COOP_REAL,parameter:: c0 = 1.d0 - 2.d0*c1
    call this%H_Split(c1*dt/2.d0,1)
    do j=1,nsteps
       call this%symp_o2step(dt, c1, c0)
       call this%symp_o2step(dt, c0, c1)
       if(j.eq.nsteps)then
          call this%symp_o2step(dt, c1, 0.d0)
       else
          call this%symp_o2step(dt, c1, c1)
       endif
    enddo
  end subroutine coop_lattice_symp4

  subroutine coop_lattice_symp6(this, dt,nsteps)
    class(coop_lattice)::this
    COOP_REAL::dt
    COOP_INT::nsteps, j
    COOP_REAL,parameter:: c1 = -1.17767998417887d0, c2 = 0.235573213359357d0, c3 = 0.784513610477560d0
    COOP_REAL,parameter:: c0 = 1.d0-2.d0*(c1+c2+c3)
    call this%H_Split(c3*dt/2.d0,1)
    do j=1,nsteps
       call this%symp_o2step(dt, c3, c2)
       call this%symp_o2step(dt, c2, c1)
       call this%symp_o2step(dt, c1, c0)
       call this%symp_o2step(dt, c0, c1)
       call this%symp_o2step(dt, c1, c2)
       call this%symp_o2step(dt, c2, c3)
       if(j.eq.nsteps)then
          call this%symp_o2step(dt, c3, 0.d0)
       else
          call this%symp_o2step(dt, c3, c3)
       endif
    enddo
  end subroutine coop_lattice_symp6

!!===================== Hamiltonian splitting ===============================
  subroutine coop_lattice_H_Split(this, dt, iterm)
    class(coop_lattice)::this
    COOP_INT::iterm
    COOP_REAL dt
    select case(iterm)
    case(1)
       call this%H_K3(dt)
    case(2)
       call this%H_K2(dt)       
    case(3)
       call this%H_K1(dt)
    case(4)
       call this%H_P(dt)
    case default
       stop "The current HLattice version only support 3-term Hamiltonian"
    end select
  end subroutine Coop_Lattice_H_Split


  function coop_lattice_K1(this) result(K1)
    class(coop_lattice)::this
    COOP_REAL::K1
    K1 = -((3.d0-this%mu)*this%pi_y)**2/48.d0/this%Mpsq
  end function coop_lattice_K1

  subroutine coop_lattice_H_K1(this, dt)
    class(coop_lattice)::this
    COOP_REAL dt
    this%y = this%y - ((3-this%mu)**2/24.d0/this%Mpsq)*this%pi_y * dt
  end subroutine coop_lattice_H_K1

  
  function coop_lattice_K2(this) result(K2)
    class(coop_lattice)::this
    COOP_REAL::K2
    COOP_INT:: LATTICE_INDS
    COOP_REAL::sqrtg
    K2 = 0.d0
    !$omp parallel do private(LATTICE_INDS) reduction(+:K2)
    LOOP_LEVEL1
    K2 = K2 &
         + sum(this%pi_f(:, imin:imax, jmin:jmax, kmin:kmax)**2)/2.d0 * exp(-sum(this%beta(1:3, ii, jj, kk))/2.d0)  &
         + sum(this%pi_beta(4:6, ii, jj, kk)**2)*(this%cs3/this%Mpsq)
    END_LOOP
    !$omp end parallel do
    K2 = K2/this%y**2
  end function coop_lattice_K2
  

  subroutine coop_lattice_H_K2(this, dt)
    class(coop_lattice)::this
    COOP_REAL dt, rdt, ebdt
    COOP_INT::LATTICE_INDS
    !!fields kinetic
    rdt = dt/this%y**2
    !$omp parallel do private(LATTICE_INDS, ebdt)
    LOOP_LEVEL1
    ebdt = exp(-sum(this%beta(1:3, ii, jj, kk))/2.d0)*rdt
    LOOP_LEVEL2
    this%f(:, i, j, k) = this%f(:, i, j, k) + ebdt * this%pi_f(:, i, j, k)
    END_LOOP
    END_LOOP
    !$omp end parallel do
    call this%wrap_f()
    
    !!background
    this%pi_y = this%pi_y + 2.d0*dt/this%y*this%K2()
    
    !!gravity kinetic off diagonal    
    !!TBD
  end subroutine coop_lattice_H_K2

  function coop_lattice_K3(this) result(K3)
    class(coop_lattice)::this
    COOP_REAL::K3
    COOP_INT::LATTICE_INDS
    !!TBD
    K3 = 0.d0
  end function coop_lattice_K3



  subroutine coop_lattice_H_K3(this, dt)
    class(coop_lattice)::this
    COOP_REAL dt
    !!gravity kinetic diagonal; TBD    
  end subroutine coop_lattice_H_K3



  function coop_lattice_P(this) result(P)
    class(coop_lattice)::this
    COOP_REAL::P
    COOP_INT::LATTICE_INDS
    P = this%G_f() + this%V_f() + this%G_g()
  end function coop_lattice_P

  subroutine coop_lattice_get_hup(this)
    class(coop_lattice)::this    
    this%hup(IND11,:,:,:) = -this%beta(IND11,:,:,:) + (this%beta(IND11,:,:,:)**2 + this%beta(IND12,:,:,:,:,:,:)**2 + this%beta(IND13,:,:,:)**2)/2.d0
    this%hup(IND22,:,:,:) = -this%beta(IND22,:,:,:) + (this%beta(IND21,:,:,:)**2 + this%beta(IND22,:,:,:)**2 + this%beta(IND23,:,:,:)**2)/2.d0
    this%hup(IND33,:,:,:) = -this%beta(IND33,:,:,:) + (this%beta(IND31,:,:,:)**2 + this%beta(IND32,:,:,:)**2 + this%beta(IND33,:,:,:)**2)/2.d0
    this%hup(IND23,:,:,:) = -this%beta(IND23,:,:,:) + (this%beta(IND21,:,:,:)*this%beta(IND31,:,:,:) + this%beta(IND22,:,:,:)*this%beta(IND32,:,:,:) + this%beta(IND23,:,:,:)*this%beta(IND33,:,:,:))/2.d0
    this%hup(IND31,:,:,:) = -this%beta(IND31,:,:,:) + (this%beta(IND31,:,:,:)*this%beta(IND11,:,:,:) + this%beta(IND32,:,:,:)*this%beta(IND12,:,:,:) + this%beta(IND33,:,:,:)*this%beta(IND13,:,:,:))/2.d0
    this%hup(IND12,:,:,:) = -this%beta(IND12,:,:,:) + (this%beta(IND11,:,:,:)*this%beta(IND21,:,:,:) + this%beta(IND12,:,:,:)*this%beta(IND22,:,:,:) + this%beta(IND13,:,:,:)*this%beta(IND23,:,:,:))/2.d0
  end subroutine coop_lattice_get_hup


  function coop_lattice_G_f(this) result(G_f)
    class(coop_lattice)::this
    COOP_REAL::G_f, eb
    COOP_INT::LATTICE_INDS
    G_f = 0.d0
    !$omp parallel do private(LATTICE_INDS, eb) reduction(+:G_f)
    LOOP_LEVEL1
    eb = exp(sum(this%beta(1:3, ii, jj, kk))/2.d0)
    LOOP_LEVEL2
    G_f = G_f  + eb * ( &
         this%gup(IND11,ii,jj,kk) * ( sum((this%f(:, i+1, j, k) - this%f(:, i, j, k))**2) + sum((this%f(:, i-1, j, k) - this%f(:, i, j, k))**2) ) &
         + this%gup(IND22,ii,jj,kk) * ( sum((this%f(:, i, j+1, k) - this%f(:, i, j, k))**2) +  sum((this%f(:, i, j-1, k) - this%f(:, i, j, k))**2) ) &
         + this%gup(IND33, ii,jj,kk) * ( sum((this%f(:, i, j, k+1) - this%f(:, i, j, k))**2) + sum((this%f(:, i, j, k-1) - this%f(:, i, j, k))**2) )) &
         + this%gup(IND23,ii,jj,kk) * sum((this%f(:, i, j+1, k) - this%f(:, i, j-1, k))*(this%f(:, i, j, k+1) - this%f(:, i, j, k-1))) &
         + this%gup(IND31,ii,jj,kk) * sum((this%f(:, i+1, j, k) - this%f(:, i-1, j, k))*(this%f(:, i, j, k+1) - this%f(:, i, j, k-1))) &
         + this%gup(IND12,ii,jj,kk) * sum((this%f(:, i, j+1, k) - this%f(:, i, j-1, k))*(this%f(:, i+1, j, k) - this%f(:, i-1, j, k))) 
    END_LOOP
    END_LOOP
    !$omp end parallel do
    G_f = G_f / 4.d0 / this%dx**2 * this%y ** (2.d0*(1.d0+this%mu)/(3.d0-this%mu))
  end function coop_lattice_G_f


  function coop_lattice_V_f(this) result(V_f)
    class(coop_lattice)::this
    COOP_REAL::V_f, eb
    COOP_INT::LATTICE_INDS
    V_f = 0.d0
    !$omp parallel do private(LATTICE_INDS, eb) reduction(+:V_f)
    LOOP_LEVEL1
    eb = exp(sum(this%beta(1:3, ii, jj, kk))/2.d0)    
    LOOP_LEVEL2
    V_f = V_f + coop_lattice_potential(this%f(:, i, j, k))*eb
    END_LOOP
    END_LOOP
    !$omp end parallel do
    V_f = V_f * this%y ** ( 2.d0*(3.d0+this%mu)/(3.d0-this%mu) )
  end function coop_lattice_V_f
  
  function coop_lattice_G_g(this) result(G_g)
    class(coop_lattice)::this
    COOP_REAL::G_g
    COOP_INT::LATTICE_INDS
    G_g = 0.d0
    !!TBD
  end function coop_lattice_G_g  

  subroutine coop_lattice_H_P(this, dt)
    class(coop_lattice)::this
    COOP_REAL dt, ebdt, yfac1, yfac2
    COOP_INT::LATTICE_INDS    
    !!potential and gradient
    yfac1 = this%y**(2.d0*(1.d0+this%mu)/(3.d0-this%mu))
    yfac2 = this%y**(2.d0*(3.d0+this%mu)/(3.d0-this%mu))
    this%pi_y = this%pi_y - (2.d0*(1.d0+this%mu)/(3.d0-this%mu)*(this%G_f()+this%G_g()) + 2.d0*(3.d0+this%mu)/(3.d0-this%mu)*this%V_f())/this%y * dt
    !$omp parallel do private(LATTICE_INDS, gup, eb)
    LOOP_LEVEL1
    ebdt = exp(sum(this%beta(1:3, ii, jj, kk))/2.d0)*dt
    call coop_lattice_beta2gup(this%beta(:, ii, jj, kk), gup)
    LOOP_LEVEL2
    this%pi_f(:, i, j, k) = this%pi_f(:, i, j, k) + ebdt * (yfac1*( &

         ) - yfac2*coop_lattice_dVdphi(this%f(:, i, j, k)))
    END_LOOP
    END_LOOP
    !$omp end parallel do

    !!metric evolution
    !!TBD
  end subroutine coop_lattice_H_P
  
  
  
  
end module coop_lattice_mod



