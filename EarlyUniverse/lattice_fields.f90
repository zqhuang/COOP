!!!lattice simulation of scalar dynamics (ignore metric perturbations)
module coop_lattice_fields_mod
  use coop_wrapper_utils
  implicit none
#include "constants.h"
#include "lattice.h"  


  private
  public:: coop_lattice_fields, coop_lattice_Mp, coop_lattice_Mpsq, coop_lattice_fields_V, coop_lattice_fields_dVdphi


  COOP_REAL, parameter::coop_lattice_Mp = 1024.d0
  COOP_REAL, parameter::coop_lattice_Mpsq = coop_lattice_Mp ** 2
!!1/(8\pi G); This just sets the program unit. 

!!************** define your parameters for potential *************
  COOP_REAL, parameter::lambda = 1.d-13 !!lambda
  COOP_REAL, parameter::g2byl = 2.d0 !!g^2/lambda
!!*****************************************************************

  
  type coop_lattice_fields
     COOP_INT::nflds = 0
     COOP_INT::n = 0
     COOP_REAL::n3 = 0.d0
     COOP_INT::n_H_terms = 3
     COOP_INT::ode_order = 6
     COOP_REAL, dimension(:,:,:,:),allocatable::f, pi
     COOP_REAL::y = 1.d0  !! y = a^((3-mu)/2)
     COOP_REAL::pi_y
     COOP_REAL::dx = 1.d0
     COOP_REAL::dk, L, kmax
     COOP_REAL::mu = 1.d0   !!d (cosmological time)  = a^mu d (program time)
     COOP_REAL::raw_KE, raw_GE, raw_PE !!not normalize
     COOP_REAL::KE, GE, PE,a, H !!normalize
     logical::expansion = .true.
   contains
     procedure::pi_y_constrained => coop_lattice_fields_pi_y_constrained
     procedure::set_pi_y => coop_lattice_fields_set_pi_y
     procedure::free => coop_lattice_fields_free
     procedure::init => coop_lattice_fields_init
     procedure::symp2 => coop_lattice_fields_symp2
     procedure::symp4 => coop_lattice_fields_symp4
     procedure::symp6 => coop_lattice_fields_symp6
     procedure::symp_o2step => coop_lattice_fields_symp_o2step
     procedure::alloc => coop_lattice_fields_alloc
     procedure::H_split => coop_lattice_fields_H_split
     procedure::H_Kinetic => coop_lattice_fields_H_Kinetic
     procedure::H_gravity => coop_lattice_fields_H_gravity
     procedure::H_Potential_and_Gradient => coop_lattice_fields_H_Potential_and_Gradient
     procedure::wrap => coop_lattice_fields_wrap
     procedure::evolve => coop_lattice_fields_evolve
     procedure::Kinetic_Energy => coop_lattice_fields_Kinetic_Energy
     procedure::Potential_Energy => coop_lattice_fields_Potential_Energy
     procedure::Gradient_Energy => coop_lattice_fields_Gradient_Energy
  end type coop_lattice_fields

contains


  !!============ editable section: define your model here ===========


#define PHI phi(1)
#define CHI phi(2)  

  function coop_lattice_fields_V(phi) result(V)
    COOP_REAL::phi(:), V
    V =  lambda * ( (1.d0/4.d0) * PHI **2 + (g2byl/2.d0) * CHI **2 ) * PHI**2
  end function coop_lattice_fields_V

  function coop_lattice_fields_dVdphi(phi) result(dVdphi)
    COOP_REAL::phi(:), dVdphi(size(phi))
    dVdphi(1) = lambda*(PHI**2 + g2byl * CHI**2) * PHI
    dVdphi(2) = (lambda*g2byl) * CHI * PHI**2
  end function coop_lattice_fields_dVdphi
  
#undef PHI
#undef CHI  
  !!===========================end of editable section ==============

  subroutine coop_lattice_fields_free(this)
    class(coop_lattice_fields)::this
    COOP_DEALLOC(this%f)
    COOP_DEALLOC(this%pi)
    this%n = 0
    this%n3 = 0.d0
    this%nflds = 0
  end subroutine coop_lattice_fields_free
  
  subroutine coop_lattice_fields_alloc(this, nflds, n)
    class(coop_lattice_fields)::this
    COOP_INT::nflds, n
    call this%free()
    this%nflds = nflds
    this%n = n
    this%n3 = dble(n)**3
    allocate(this%f(nflds, -1:n, -1:n, -1:n))
    allocate(this%pi(nflds, 0:n-1, 0:n-1, 0:n-1))
  end subroutine coop_lattice_fields_alloc

  subroutine coop_lattice_fields_wrap(this)
    class(coop_lattice_fields)::this    
    this%f(:, -1, :, :) = this%f(:, this%n - 1, :, :)
    this%f(:, :, -1, :) = this%f(:, :, this%n - 1, :)
    this%f(:, :, :, -1) = this%f(:, :, :, this%n - 1)
    this%f(:, this%n, :, :) = this%f(:, 0, :, :)
    this%f(:, :, this%n, :) = this%f(:, :, 0, :)
    this%f(:, :, :, this%n) = this%f(:, :, :, 0)
  end subroutine coop_lattice_fields_wrap



  subroutine coop_lattice_fields_evolve(this, dt, nsteps)
    class(coop_lattice_fields)::this
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
  end subroutine coop_lattice_fields_evolve

  subroutine coop_lattice_fields_symp_o2step(this, dt,c1,c2)
    class(coop_lattice_fields)::this
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
  end subroutine coop_lattice_fields_symp_o2step

  subroutine coop_lattice_fields_symp2(this, dt,nsteps)
    class(coop_lattice_fields)::this
    COOP_REAL::dt
    COOP_INT::nsteps, j
    call this%H_split(dt/2.d0, 1)
    do j=1,nsteps-1
       call this%symp_o2step(dt, 1.d0, 1.d0)
    enddo
    call this%symp_o2step(dt, 1.d0, 0.d0)
  end subroutine coop_lattice_fields_symp2

  subroutine coop_lattice_fields_symp4(this, dt,nsteps)
    class(coop_lattice_fields)::this
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
  end subroutine coop_lattice_fields_symp4

  subroutine coop_lattice_fields_symp6(this, dt,nsteps)
    class(coop_lattice_fields)::this
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
  end subroutine coop_lattice_fields_symp6

!!===================== Hamiltonian splitting ===============================
  subroutine coop_lattice_fields_H_Split(this, dt, iterm)
    class(coop_lattice_fields)::this
    COOP_INT::iterm
    COOP_REAL dt
    select case(iterm)
    case(1)
       call this%H_kinetic(dt)
    case(2)
       call this%H_gravity(dt)
    case(3)
       call this%H_Potential_and_Gradient(dt)
    case default
       stop "HLattice without metric perturbations only support 3-term Hamiltonian"
    end select
  end subroutine Coop_lattice_fields_H_Split


  function coop_lattice_fields_gradient_energy(this, normalize) result(ge)
    class(coop_lattice_fields)::this
    COOP_REAL::ge
    logical,optional::normalize
    COOP_INT::i, j, k
    ge = 0.d0
    !$omp parallel do private(i, j, k) reduction(+:ge)
    DO_LOOP
    ge = ge + sum((this%f(:,i, j, k) - this%f(:,i+1, j, k))**2 &
         + (this%f(:, i, j, k) - this%f(:, i, j+1, k))**2 &
         + (this%f(:, i, j, k) - this%f(:, i, j, k+1))**2)
    END_LOOP
    ge = ge/2.d0/this%dx**2/this%n3
    if(present(normalize))then
       if(normalize)then
          ge = ge / this%y**(4.d0/(3.d0-this%mu))
       endif
    else
       ge = ge / this%y**(4.d0/(3.d0-this%mu))
    endif
  end function coop_lattice_fields_gradient_energy

  function coop_lattice_fields_potential_energy(this) result(V)
    class(coop_lattice_fields)::this
    COOP_REAL::V
    COOP_INT::i, j, k
    V = 0.d0
    !$omp parallel do private(i, j, k) reduction(+:V)
    DO_LOOP
    V = V + coop_lattice_fields_V(this%f(:,i, j, k))
    END_LOOP
    !$omp end parallel do
    V = V/this%n3
  end function coop_lattice_fields_potential_energy


  function coop_lattice_fields_kinetic_energy(this, normalize) result(Kin)
    class(coop_lattice_fields)::this
    COOP_REAL::Kin
    logical, optional::normalize
    COOP_INT::i, j, k
    Kin = sum(this%pi**2)/2.d0/this%n3
    if(present(normalize))then
       if(normalize)then
          Kin = Kin/this%y**(12.d0/(3.d0-this%mu))
       endif
    else
       Kin = Kin/this%y**(12.d0/(3.d0-this%mu))
    endif
  end function coop_lattice_fields_kinetic_energy



  subroutine coop_lattice_fields_H_gravity(this, dt)
    class(coop_lattice_fields)::this
    COOP_REAL dt
    if(this%expansion)this%y = this%y + GR_COEF * this%pi_y * dt
  end subroutine coop_lattice_fields_H_gravity


  subroutine coop_lattice_fields_H_Kinetic(this, dt)
    class(coop_lattice_fields)::this
    COOP_REAL dt
    COOP_INT::LATTICE_INDS
    F_ALL = F_ALL + this%pi*(dt * this%y**K_INDEX)
    call this%wrap()
    
    if(this%expansion) this%pi_y = this%pi_y - (K_INDEX * this%y**(K_INDEX - 1.d0)*dt ) *this%Kinetic_Energy(normalize = .false.)

  end subroutine coop_lattice_fields_H_Kinetic


  subroutine coop_lattice_fields_H_Potential_and_Gradient(this, dt)
    class(coop_lattice_fields)::this
    COOP_REAL dt, ebdt, yfac1, yfac2
    COOP_INT::i, j, k
    !!potential and gradient
    yfac1 = this%y**G_INDEX/this%dx**2*dt
    yfac2 = this%y**V_INDEX*dt
    if(this%expansion) this%pi_y = this%pi_y - (G_INDEX * this%y**(G_INDEX-1.d0) * dt) * this%Gradient_Energy(normalize = .false.) - (V_INDEX * this%y**(V_INDEX-1.d0) * dt) * this%Potential_Energy()
    !$omp parallel do private(i, j, k)
    DO_LOOP
    this%pi(:, i, j, k) = this%pi(:, i, j, k) + yfac1 * LAP_F(i, j, k) - yfac2*coop_lattice_fields_dVdphi(this%f(:, i, j, k))
    END_LOOP
    !$omp end parallel do
  end subroutine coop_lattice_fields_H_Potential_and_Gradient
  
  
  function coop_lattice_fields_pi_y_constrained(this) result(piy)
    class(coop_lattice_fields)::this
    COOP_REAL::piy
    this%raw_KE = this%kinetic_energy(normalize = .false.)
    this%KE = this%raw_KE*this%y**(-12.d0/(3.d0-this%mu))
    this%raw_PE = this%potential_energy()
    this%raw_GE = this%gradient_energy(normalize = .false.)
    this%a = this%y ** A_INDEX
    this%KE = this%raw_KE / this%a**6
    this%GE = this%raw_GE /this%a**2
    this%PE = this%raw_PE
    this%H = sqrt((-this%pi_y**2/2.d0*GR_COEF)/this%y**(2.d0*(3.d0+this%mu)/(3.d0-this%mu))/3.d0/coop_lattice_Mpsq)
    piy = - sqrt(-(this%raw_KE*this%y**K_INDEX + this%raw_GE*this%y**G_INDEX+this%raw_PE*this%y**V_INDEX)*2.d0/GR_COEF)
  end function coop_lattice_fields_pi_y_constrained


  subroutine coop_lattice_fields_set_pi_y(this)
    class(coop_lattice_fields)::this
    this%pi_y = this%pi_y_constrained()
  end subroutine coop_lattice_fields_set_pi_y

  subroutine coop_lattice_fields_init(this, mu, n, LH,  phi, pi, phi_sigma2)
    class(coop_lattice_fields)::this
    COOP_REAL::LH, phi(:), pi(:), phi_sigma2(:)
    COOP_REAL,optional::mu
    COOP_INT::n
    COOP_COMPLEX fk(0:n/2,0:n-1, 0:n-1)
    COOP_REAL::ftmp(0:n-1, 0:n-1, 0:n-1), norm, phi_sigma, k2
    COOP_INT::i, j, k, fld
    call this%alloc(nflds = size(phi), n = n)
    if(present(mu))this%mu = mu
    this%L = LH/sqrt((coop_lattice_fields_V(phi) + sum(pi**2)/2.d0)/3.d0/coop_lattice_Mpsq)
    this%dx = this%L / this%n
    this%dk = coop_2pi / this%L
    this%kmax = this%n/2 * this%dk
    norm = 1.d0/sqrt(4.d0*coop_pi)
    do fld = 1, this%nflds
       this%pi(fld, :,:,:) = pi(fld)
       phi_sigma = sqrt(phi_sigma2(fld))*norm
       do  k = 0, this%n-1; do j = 0, this%n-1; do i = 0, this%n/2
          k2 = (min(dble(j), dble(this%n-j))**2  + min(dble(k), dble(this%n-k))**2 + dble(i)**2)
          if(k2 .gt. 0.d0 .and. k2 .lt. (this%n*0.49d0)**2)then
             fk(i, j, k) = phi_sigma/k2**0.75d0*norm*coop_random_complex_Gaussian()
          else
             fk(i, j, k) = 0.d0
          endif
       enddo;enddo;enddo
       call fft_3d_backward(n, n, n, fk, ftmp)
       this%f(fld, 0:this%n-1, 0:this%n-1, 0:this%n-1) = ftmp + phi(fld)
    enddo
    call this%wrap()
    this%y = 1.d0
    call this%set_pi_y()
  end subroutine coop_lattice_fields_init


  
end module coop_lattice_fields_mod



