!! HLattice: lattice simulation of scalar dynamics (ignore metric perturbations)
!! The model is defined in include/lattice_field_model.h (For your own model, replace the file and recompile HLattice)

module coop_lattice_fields_mod
  use coop_wrapper_utils
  implicit none
#include "constants.h"
#include "lattice.h"  

  private
  public:: coop_lattice_fields, coop_lattice_Mp, coop_lattice_Mpsq, coop_lattice_fields_V, coop_lattice_fields_dVdphi, coop_lattice_background_eqs,  coop_lattice_background_epsilon,  coop_lattice_background_rhotot, coop_lattice_initial_power, coop_inflation_background, coop_infbg
  

  !!M_p^2  = 1/(8\pi G); This just sets the program unit and in principle can be arbitrary.
  COOP_REAL, parameter::coop_lattice_Mp = 1024.d0
  COOP_REAL, parameter::coop_lattice_Mpsq = coop_lattice_Mp ** 2
  COOP_REAL,parameter::Nmatter = 0. !!number of efolds of matter dominant at the end of inflation; this will affect the conversion between kMpc and aH
  !!************** define your parameters for potential *************
  !  COOP_REAL,parameter::mphi = 1.21d-5*coop_lattice_Mp !!M
  COOP_REAL,parameter::lambda_phi = 1.d-6
  COOP_INT,parameter::n_phi = 6
  COOP_REAL,parameter::alpha = 1.d-5 !!amplitude of bump
  COOP_REAL,parameter::phipk = 5.37d0*coop_lattice_Mp !!emm need to adjust later 
  COOP_REAL,parameter::mu = 0.005*coop_lattice_Mp  !!length of bump
  COOP_REAL,parameter::sigma = 3.d-8*coop_lattice_Mp  !!width of wall
  !!*****************************************************************


  type coop_inflation_background
     COOP_INT::nflds = 0
     COOP_INT::nsteps = 0
     COOP_REAL,dimension(:),allocatable::lnH, eps, lna, lnH2, eps2
     COOP_REAL,dimension(:,:),allocatable::f, fd, fdd, f2, fd2, fdd2, m2mat, m2mat2
     COOP_REAL::lnHend, nefolds
   contains
     procedure::free => coop_inflation_background_free
     procedure::alloc => coop_inflation_background_alloc
     procedure::setup =>coop_inflation_background_setup
     procedure::Hubble => coop_inflation_background_Hubble
     procedure::epsilon => coop_inflation_background_epsilon
     procedure::dlnepsdlna => coop_inflation_background_dlnepsdlna     
     procedure::fields => coop_inflation_background_fields
     procedure::dot_fields => coop_inflation_background_dot_fields
     procedure::ddot_fields => coop_inflation_background_ddot_fields     
     procedure::field => coop_inflation_background_field
     procedure::dot_field => coop_inflation_background_dot_field
     procedure::ddot_field => coop_inflation_background_ddot_field
     procedure::msqs => coop_inflation_background_msqs
     procedure::msq => coop_inflation_background_msq
     procedure::lnaH => coop_inflation_background_lnaH
     procedure::lnkMpc => coop_inflation_background_lnkMpc     
     procedure::lna_of_lnk => coop_inflation_background_lna_of_lnk     
  end type coop_inflation_background

  type coop_lattice_initial_power
     COOP_INT::nk = 0
     COOP_INT::nflds = 0
     logical::is_diagonal = .true.
     COOP_REAL,dimension(:),allocatable::lnk
     COOP_COMPLEX,dimension(:, :),allocatable::fdbyf
     COOP_REAL,dimension(:,:,:),allocatable::f_cov
   contains
     procedure::free => coop_lattice_initial_power_free
     procedure::alloc => coop_lattice_initial_power_alloc     
     procedure::init => coop_lattice_initial_power_initialize
  end type coop_lattice_initial_power


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
     COOP_REAL::mu = 1.d0
     !!d (cosmological time)  = a^mu d (program time);
     COOP_REAL::raw_KE, raw_GE, raw_PE !!not normalize
     COOP_REAL::KE, GE, PE,a, H !!normalize
     logical::expansion = .true.
   contains
     procedure::scale_factor => coop_lattice_fields_scale_factor
     procedure::pi_y_constrained => coop_lattice_fields_pi_y_constrained
     procedure::set_pi_y => coop_lattice_fields_set_pi_y
     procedure::set_energies => coop_lattice_fields_set_energies     
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


  type(coop_inflation_background)::coop_infbg

contains


  !!============ include model file =====================
#include "lattice_field_model.h"
  !!------------------------------------

  
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
    call this%set_energies()
    piy = - sqrt(-(this%raw_KE*this%y**K_INDEX + this%raw_GE*this%y**G_INDEX+this%raw_PE*this%y**V_INDEX)*2.d0/GR_COEF)
  end function coop_lattice_fields_pi_y_constrained


  subroutine coop_lattice_fields_set_pi_y(this)
    class(coop_lattice_fields)::this
    this%pi_y = this%pi_y_constrained()
  end subroutine coop_lattice_fields_set_pi_y


  subroutine coop_lattice_fields_set_energies(this)
    class(coop_lattice_fields)::this
    this%raw_KE = this%kinetic_energy(normalize = .false.)
    this%KE = this%raw_KE*this%y**(-12.d0/(3.d0-this%mu))
    this%raw_PE = this%potential_energy()
    this%raw_GE = this%gradient_energy(normalize = .false.)
    this%a = this%y ** A_INDEX
    this%KE = this%raw_KE / this%a**6
    this%GE = this%raw_GE /this%a**2
    this%PE = this%raw_PE
    this%H = sqrt((-this%pi_y**2/2.d0*GR_COEF)/this%y**(2.d0*(3.d0+this%mu)/(3.d0-this%mu))/3.d0/coop_lattice_Mpsq)
  end subroutine coop_lattice_fields_set_energies

  function coop_lattice_fields_scale_factor(this) result(a)
    class(coop_lattice_fields)::this
    COOP_REAL::a
    a = this%y**(2.d0/(3.d0-this%mu))
  end function coop_lattice_fields_scale_factor

  !!set initial perturbatioins
  !!n: the box resolultion (n^3 grids)
  !!LH: the comoving boxsize in unit of (aH)^{-1}
  !!phi: initial background fields
  !!pi: initial background d phi/ dt
  !!
  !!use_conformal_time (default false): if true, use conformal time;  if false, use physical time; 
  subroutine coop_lattice_fields_init(this, n, LH,  phi, pi, inipower, use_conformal_time)
    class(coop_lattice_fields)::this
    class(coop_lattice_initial_power)::inipower
    COOP_REAL :: LH, phi(:), pi(:)
    COOP_INT::n
    COOP_REAL::norm,  k2
    COOP_COMPLEX,dimension(:,:,:,:), allocatable::fk, fdk
    COOP_INT::i, j, k, fld, ik
    logical,optional::use_conformal_time
    COOP_REAL::lnk, lndk, rk, Hubble

    call this%alloc(nflds = size(phi), n = n)
    if(present(use_conformal_time))then
       if(use_conformal_time)then
          this%mu = 1.d0
       else
          this%mu = 0.d0
       endif
    else
       this%mu = 0.d0
    endif

    if(  size(pi) .ne. this%nflds)call coop_return_error("lattice_fields_init", "array sizes do not match", "stop")
    Hubble = sqrt((coop_lattice_fields_V(phi) + sum(pi**2)/2.d0)/3.d0/coop_lattice_Mpsq)
    this%L = LH/Hubble
    this%dx = this%L / this%n
    this%dk = coop_2pi / this%L
    lndk = log(this%dk)
    this%kmax = this%n/2 * this%dk
    norm = 1.d0/sqrt(4.d0*coop_pi)
    if(inipower%is_diagonal)then
       allocate(fk(0:this%n/2, 0:this%n-1, 0:this%n-1, 1), fdk(0:this%n/2, 0:this%n-1, 0:this%n-1, 1))
       do fld = 1, this%nflds
          do  k = 0, this%n-1; do j = 0, this%n-1; do i = 0, this%n/2
             k2 = (min(dble(j), dble(this%n-j))**2  + min(dble(k), dble(this%n-k))**2 + dble(i)**2)
             lnk = log(k2)/2.d0 + lndk
             if(k2 .gt. 0.d0 .and. k2 .lt. (this%n*0.49999d0)**2)then
                ik = coop_left_index(n, inipower%lnk, lnk)
                if(ik .gt. 0 .and. ik .lt. inipower%nk)then
                   rk = (inipower%lnk(ik+1)-lnk)/(inipower%lnk(ik+1)-inipower%lnk(ik))
                   fk(i, j, k, 1) = norm*sqrt(inipower%f_cov(fld, fld, ik)*rk + inipower%f_cov(fld, fld, ik+1)*(1.d0-rk)) / k2**0.75d0 *  coop_random_complex_Gaussian()
                   fdk(i, j, k, 1) =  fk(i, j, k, 1) * (inipower%fdbyf(fld, ik)*rk + inipower%fdbyf(fld, ik+1)*(1.d0-rk))
                else
                   fk(i, j, k,1) = 0.d0
                   fdk(i, j, k,1) = 0.d0                   
                endif
             else
                fk(i, j, k, 1) = 0.d0
                fdk(i, j, k,1) = 0.d0                                   
             endif
          enddo;enddo;enddo
          fk(0,0,0,1) = phi(fld)
          fdk(0,0,0,1) = pi(fld)
          call fft_3d_backward(n, n, n, fk(:,:,:,1), this%f(fld, 0:this%n-1, 0:this%n-1, 0:this%n-1))
          call fft_3d_backward(n, n, n, fdk(:,:,:,1), this%pi(fld, 0:this%n-1, 0:this%n-1, 0:this%n-1))          
       enddo
       deallocate(fk, fdk)
    else
       stop "Non-diagonal initial conditions are not yet implemented."
    endif
    call this%wrap()
    this%y = 1.d0
    call this%set_pi_y()
  end subroutine coop_lattice_fields_init

  !!================ utilities for background evolution ====================  
  ! n = 2 * nflds + 2
  !y: phi, dphi/dt, lna, H
  !yp: dy /d t
  subroutine coop_lattice_background_eqs(n, t, y, yp)
    COOP_INT::n
    COOP_REAL::t, y(n), yp(n)
    COOP_INT::m
    m = n/2-1
    yp(1:m) = y(m+1:2*m)
    yp(m+1:2*m) = -3.d0*y(m+1:2*m)*y(n) - coop_lattice_fields_dVdphi(y(1:m))
    yp(n-1) = y(n)
    yp(n) = (-0.5d0/coop_lattice_Mpsq) * sum(y(m+1:2*m)**2)
  end subroutine coop_lattice_background_eqs

  !!epsilon = -\dot H / H^2
  function coop_lattice_background_epsilon(nflds, y) result(eps)
    COOP_INT::nflds
    COOP_REAL::eps, y(2*nflds+2)
    eps = (0.5d0/coop_lattice_Mpsq) * sum(y(nflds+1:2*nflds)**2)/y(2*nflds+2)**2
  end function coop_lattice_background_epsilon


  !!total density
  function coop_lattice_background_rhotot(nflds, y) result(rho)
    COOP_INT::nflds
    COOP_REAL::rho, y(2*nflds+2)
    rho =  sum(y(nflds+1:2*nflds)**2)/2.d0 + coop_lattice_fields_V(y(1:nflds))
  end function coop_lattice_background_rhotot

  subroutine coop_inflation_background_free(this)
    class(coop_inflation_background)::this
    if(this%nsteps .gt. 0 .and. this%nflds .gt. 0)then
       deallocate(this%lna, this%lnH, this%eps, this%lnH2, this%eps2, this%f, this%fd, this%fdd, this%f2, this%fd2, this%fdd2, this%m2mat, this%m2mat2)
     this%nsteps = 0
     this%nflds = 0
    endif
  end subroutine coop_inflation_background_free


  subroutine coop_inflation_background_alloc(this, nflds, nsteps)
    class(coop_inflation_background)::this
    COOP_INT::nflds, nsteps
    if(this%nflds .ne. nflds .or. this%nsteps .ne. this%nsteps)then
       call this%free()
       allocate(this%lna(nsteps), this%lnH(nsteps), this%eps(nsteps), this%lnH2(nsteps), this%eps2(nsteps), this%f(nsteps, nflds), this%fd(nsteps,nflds), this%fdd(nsteps, nflds), this%f2(nsteps, nflds), this%fd2(nsteps, nflds), this%fdd2(nsteps, nflds), this%m2mat(nsteps, nflds*(nflds+1)/2), this%m2mat2(nsteps, nflds*(nflds+1)/2))
       this%nflds = nflds
       this%nsteps = nsteps       
    endif
  end subroutine coop_inflation_background_alloc

  subroutine coop_inflation_background_setup(this, nflds, epsilon_end, f_ini, fd_ini)
    class(coop_inflation_background)::this
    COOP_INT::nflds
    COOP_REAL::f_ini(nflds)
    COOP_REAL,optional::fd_ini(nflds)
    COOP_REAL::epsilon_end
    type(coop_ode)::bg
    !!other variables
    COOP_REAL::Hini, dotf_ini(nflds)
    COOP_REAL::y(nflds*2+2), dt, Vpp(nflds, nflds)
    COOP_INT::i, j, istep, k
    type(coop_list_double)::lnH, lna, eps
    type(coop_list_doublearr)::f, fd, fdd
    logical::eps_increase = .false.

    !!===============initialize ODE solver============================
    call bg%init(n=2*nflds+2, method=COOP_ODE_DVERK, tol=1.d-8)
    !!=================set inital conditions =========================================
    if(present(fd_ini))then
       dotf_ini = fd_ini
       Hini = sqrt(coop_lattice_fields_V(f_ini) + sum(dotf_ini**2)/2.d0)/(coop_sqrt3*coop_lattice_Mp)
    else
       !!slow-roll approximation
       Hini = sqrt(coop_lattice_fields_V(f_ini))/(coop_sqrt3*coop_lattice_Mp)
       dotf_ini  = coop_lattice_fields_dVdphi(f_ini)/(-3.d0*Hini)
       if(sum(dotf_ini**2) .gt. 0.5*coop_lattice_fields_V(f_ini) )then
          write(*,*) "Warning: slow-roll condition is not satisfied"
          write(*,*) "--------------------------------------------------" 
          write(*,*) sum(dotf_ini**2)/2.d0, coop_lattice_fields_V(f_ini)
          write(*,*) "--------------------------------------------------"     
       else
          !!iterate to get more accurate values
          Hini = sqrt(coop_lattice_fields_V(f_ini) + sum(dotf_ini**2)/2.d0)/(coop_sqrt3*coop_lattice_Mp)
          dotf_ini  = coop_lattice_fields_dVdphi(f_ini)/(-3.d0*Hini)
          Hini = sqrt(coop_lattice_fields_V(f_ini) + sum(dotf_ini**2)/2.d0)/(coop_sqrt3*coop_lattice_Mp)
       endif
    endif
    y(1:nflds) = f_ini
    y(nflds+1:2*nflds) = dotf_ini
    y(2*nflds+1) = 0.d0  !!initial ln a = 0
    y(2*nflds+2) = Hini
    call bg%set_initial_conditions(xini = 0.d0, yini = y)
    !!======================evolution=================================================
#define LNA bg%y(2*nflds+1)
#define HUBBLE bg%y(2*nflds+2)
#define FIELDS bg%y(1:nflds)
#define DOT_FIELDS    bg%y(nflds+1:2*nflds)
    dt = 0.1/Hini
    call lna%push(LNA)
    call lnH%push(log(HUBBLE))
    call eps%push(coop_lattice_background_epsilon(nflds, bg%y))
    call f%push(FIELDS)
    call fd%push(DOT_FIELDS)
    call coop_lattice_background_eqs(bg%n, bg%x, bg%y, y)
    call fdd%push(y(nflds+1:2*nflds))
    
    if(eps%element(1) .gt. epsilon_end) call coop_return_error("coop_inflation_background_setup",  "Initial epsilon overflow.", "stop")
    do while(eps%element(eps%n) .lt. epsilon_end*0.9999d0)
       if(eps_increase .and. eps%element(eps%n) .gt.  epsilon_end*0.5d0 )then
          dt = max(min(0.3 * dt * (epsilon_end-eps%element(eps%n))/(eps%element(eps%n) - eps%element(eps%n-1)), 0.1d0/HUBBLE), 1.d-6/HUBBLE)
       else
          dt = 0.1/HUBBLE
       endif
       call bg%evolve(coop_lattice_background_eqs, bg%x+dt)
       call lna%push(LNA)
       call lnH%push(log(HUBBLE))
       call eps%push(coop_lattice_background_epsilon(nflds, bg%y))
       call f%push(FIELDS)
       call fd%push(DOT_FIELDS)
       call coop_lattice_background_eqs(bg%n, bg%x, bg%y, y)
       call fdd%push(y(nflds+1:2*nflds))
       eps_increase = (eps%element(eps%n) .gt. eps%element(eps%n-1) + 1.d-8)

       if(abs(coop_lattice_background_rhotot(nflds, bg%y) /(3.d0*coop_lattice_Mpsq*HUBBLE**2) - 1.d0) .gt. 1.d-4) stop "energy conservation failed: check if you have a typo in include/lattice_fields_mode.h"

    enddo
    dt = dt * (epsilon_end - eps%element(eps%n))/(eps%element(eps%n) - eps%element(eps%n-1))
    call bg%evolve(coop_lattice_background_eqs, bg%x+dt)
    call lna%push(LNA)
    call lnH%push(log(HUBBLE))
    call eps%push(coop_lattice_background_epsilon(nflds, bg%y))
    call f%push(FIELDS)
    call fd%push(DOT_FIELDS)
    call coop_lattice_background_eqs(bg%n, bg%x, bg%y, y)
    call fdd%push(y(nflds+1:2*nflds))
    
#undef LNA
#undef HUBBLE
#undef FIELDS
#undef DOT_FIELDS
    call this%alloc(nflds = nflds, nsteps = lna%n)
    do istep=1, lna%n
       call lna%get_element(istep, this%lna(istep))
       call lnH%get_element(istep, this%lnH(istep))
       call eps%get_element(istep, this%eps(istep))
       this%f(istep,:) = f%element(istep)
       this%fd(istep,:) = fd%element(istep)
       this%fdd(istep, :) = fdd%element(istep)
       Vpp =  coop_lattice_fields_d2Vdphi2(this%f(istep,:))
       do i=1, nflds
          do j=1, i
             this%m2mat(istep, COOP_MATSYM_INDEX(nflds, i, j))  = Vpp(i, j) - (this%fd(istep, i)*this%fd(istep,j)*(3.d0+this%eps(istep)) + (this%fd(istep,i)*this%fdd(istep,j) + this%fd(istep,j)*this%fdd(istep,i))/exp(this%lnH(istep)))/coop_lattice_Mpsq
             Vpp(j, i) = Vpp(i, j)
          enddo
       enddo
    enddo
    this%lna = this%lna - this%lna(this%nsteps)
    call coop_spline(this%nsteps, this%lna, this%eps, this%eps2)
    call coop_spline(this%nsteps, this%lna, this%lnH, this%lnH2)
    do i=1, nflds
       call coop_spline(this%nsteps, this%lna, this%f(:,i), this%f2(:,i))
       call coop_spline(this%nsteps, this%lna, this%fd(:,i), this%fd2(:,i))
       call coop_spline(this%nsteps, this%lna, this%fdd(:,i), this%fdd2(:,i))  
    enddo
    do i=1, nflds*(nflds+1)/2
       call coop_spline(this%nsteps, this%lna, this%m2mat(:, i), this%m2mat2(:, i))
    enddo
    this%lnHend = this%lnH(this%nsteps)
    this%nefolds = - this%lna(1)
    call bg%free()
    call lna%free()
    call lnH%free()
    call eps%free()
    call f%free()
    call fd%free()
    call fdd%free()
  end subroutine coop_inflation_background_setup

  function coop_inflation_background_Hubble(this, lna) result(H)
    class(coop_inflation_background)::this
    COOP_REAL::lna, H
    call coop_splint(this%nsteps, this%lna, this%lnH, this%lnH2, lna, H)
    H = exp(H)
  end function coop_inflation_background_Hubble

  function coop_inflation_background_lnaH(this, lna) result(lnaH)
    class(coop_inflation_background)::this
    COOP_REAL::lna, lnaH
    call coop_splint(this%nsteps, this%lna, this%lnH, this%lnH2, lna, lnaH)
    lnaH = lnaH + lna
  end function coop_inflation_background_lnaH

  function coop_inflation_background_lna_of_lnk(this, lnk) result(lna)
    class(coop_inflation_background)::this
    COOP_REAL::lna, lnk, bottom, top
    bottom = this%lna(1)
    top = this%lna(this%nsteps)
    do while(top - bottom .gt. 1.d-6)
       lna = (bottom + top)/2.d0
       if(this%lnaH(lna) .lt. lnk)then
          bottom = lna
       else
          top = lna
       endif
    enddo
  end function coop_inflation_background_lna_of_lnk


  function coop_inflation_background_lnkMpc(this) result(lnkMpc)
    class(coop_inflation_background)::this    
    COOP_REAL::lnkMpc
    lnkMpc = log(3.5d-26*sqrt(coop_lattice_Mp))+0.5*this%lnHend + 0.25*Nmatter
  end function coop_inflation_background_lnkMpc
  
  function coop_inflation_background_epsilon(this, lna) result(eps)
    class(coop_inflation_background)::this
    COOP_REAL::lna, eps
    call coop_splint(this%nsteps, this%lna, this%eps, this%eps2, lna, eps)
  end function coop_inflation_background_epsilon


  function coop_inflation_background_fields(this, lna) result(f)
    class(coop_inflation_background)::this
    COOP_REAL::lna, f(this%nflds)
    COOP_INT::i
    do i=1, this%nflds
       call coop_splint(this%nsteps, this%lna, this%f(:, i), this%f2(:, i), lna, f(i))
    enddo
  end function coop_inflation_background_fields


  function coop_inflation_background_dot_fields(this, lna) result(fd)
    class(coop_inflation_background)::this
    COOP_REAL::lna, fd(this%nflds)
    COOP_INT::i
    do i=1, this%nflds
       call coop_splint(this%nsteps, this%lna, this%fd(:, i), this%fd2(:, i), lna, fd(i))
    enddo
  end function coop_inflation_background_dot_fields


  function coop_inflation_background_ddot_fields(this, lna) result(fdd)
    class(coop_inflation_background)::this
    COOP_REAL::lna, fdd(this%nflds)
    COOP_INT::i
    do i=1, this%nflds
       call coop_splint(this%nsteps, this%lna, this%fdd(:, i), this%fdd2(:, i), lna, fdd(i))
    enddo
  end function coop_inflation_background_ddot_fields

  function coop_inflation_background_dlnepsdlna(this, lna) result(dlnepsdlna)
    class(coop_inflation_background)::this
    COOP_REAL::lna, dlnepsdlna, eps, df(this%nflds), ddf(this%nflds), H
    df = this%dot_fields(lna)
    ddf = this%ddot_fields(lna)
    H = this%Hubble(lna)
    eps = sum(df**2)/H**2/(2.d0*coop_lattice_Mpsq)
    dlnepsdlna = sum(df*ddf)/eps/H**3/coop_lattice_Mpsq + 2.d0*eps
  end function coop_inflation_background_dlnepsdlna

  function coop_inflation_background_msqs(this, lna) result(msqs)
    class(coop_inflation_background)::this
    COOP_REAL::lna, msqs(this%nflds, this%nflds)
    COOP_INT::i, j, k
    do i=1, this%nflds
       do j=1, i
          k = COOP_MATSYM_INDEX(this%nflds, i, j)
          call coop_splint(this%nsteps, this%lna, this%m2mat(:, k), this%m2mat2(:, k), lna, msqs(i,j))
          msqs(j, i) = msqs(i, j)
       enddo
    enddo
  end function coop_inflation_background_msqs


    function coop_inflation_background_field(this, lna, i) result(f)
    class(coop_inflation_background)::this
    COOP_REAL::lna, f
    COOP_INT::i
    call coop_splint(this%nsteps, this%lna, this%f(:, i), this%f2(:, i), lna, f)
  end function coop_inflation_background_field


  function coop_inflation_background_dot_field(this, lna, i) result(fd)
    class(coop_inflation_background)::this
    COOP_REAL::lna, fd
    COOP_INT::i
    call coop_splint(this%nsteps, this%lna, this%fd(:, i), this%fd2(:, i), lna, fd)
  end function coop_inflation_background_dot_field


  function coop_inflation_background_ddot_field(this, lna, i) result(fdd)
    class(coop_inflation_background)::this
    COOP_REAL::lna, fdd
    COOP_INT::i
    call coop_splint(this%nsteps, this%lna, this%fd(:, i), this%fd2(:, i), lna, fdd)
  end function coop_inflation_background_ddot_field
  
  

  function coop_inflation_background_msq(this, lna, i, j) result(msq)
    class(coop_inflation_background)::this
    COOP_REAL::lna, msq
    COOP_INT::i, j
    call coop_splint(this%nsteps, this%lna, this%m2mat(:, COOP_MATSYM_INDEX(this%nflds, i, j)), this%m2mat2(:, COOP_MATSYM_INDEX(this%nflds, i, j)), lna, msq)
  end function coop_inflation_background_msq
  
  

  !!================ utilities for perturbation evolution ====================  

  subroutine coop_lattice_initial_power_free(this)
    class(coop_lattice_initial_power)::this
    COOP_DEALLOC(this%lnk)
    COOP_DEALLOC(this%f_cov)
    COOP_DEALLOC(this%fdbyf)
    this%nk = 0
    this%nflds = 0
  end subroutine coop_lattice_initial_power_free


  subroutine coop_lattice_initial_power_alloc(this, nk, nflds)
    class(coop_lattice_initial_power)::this
    COOP_INT::nk, nflds
    call this%free()
    this%nk = nk
    this%nflds = nflds
    allocate(this%lnk(nk), this%f_cov(nflds, nflds, nk), this%fdbyf(nflds, nk))
  end subroutine coop_lattice_initial_power_alloc
  

  subroutine coop_lattice_initial_power_initialize(this, lnkmin, lnkmax, is_diagonal)
    COOP_INT,parameter::n_trials = 100
    COOP_REAL,parameter::shift = log(1000.d0) !!shift so many efolds ahead to initialize subhorizon modes
    class(coop_lattice_initial_power)::this
    logical, optional::is_diagonal
    COOP_REAL,optional::lnkmin, lnkmax
    COOP_REAL::omega
    type(coop_ode)::pert_r, pert_i
    COOP_REAL::ini_r(2*coop_infbg%nflds), ini_i(2*coop_infbg%nflds), lna, msq, kbya, kbyasq, maxmsq, amp, theta
    COOP_INT::ik, i, j, n_seeds, fld, i_start
    type(coop_arguments)::args
    if(coop_infbg%nflds .eq. 0 .or. coop_infbg%nsteps .eq. 0) call coop_return_error("coop_lattice_initial_power_initialize", "You need to set up inflation background before calculating perturbations", "stop")
    if(present(lnkmin) .and. present(lnkmax))then
       if(coop_infbg%lna(1)+coop_infbg%lnH(1) .ge. lnkmin - shift*0.99d0) call coop_return_error("coop_lattice_initial_power_initialize", "kmin is too small, you need to calculate more efolds for background", "stop")
       call coop_set_uniform(this%nk, this%lnk, lnkmin, lnkmax)
    endif
    if(present(is_diagonal))then
       this%is_diagonal = is_diagonal
    else
       this%is_diagonal = .true.
    endif
    if(this%is_diagonal)then
       this%f_cov = 0.
       this%fdbyf = 0.
       call pert_r%init(n=2, method=COOP_ODE_DVERK, tol = 1.d-8)
       call pert_i%init(n=2, method=COOP_ODE_DVERK, tol = 1.d-8)
       do ik = 1, this%nk
          i = 1
          do while(coop_infbg%lna(i)+coop_infbg%lnH(i) .lt. this%lnk(ik) - shift .and. i.lt. coop_infbg%nsteps)
             i = i+1
          enddo
          i_start = i - 1
          lna = coop_infbg%lna(i_start) 
          kbya = exp(this%lnk(ik)-lna)
          kbyasq = kbya ** 2
          do fld = 1, coop_infbg%nflds
             maxmsq = kbyasq
             msq = coop_infbg%msq(lna, fld, fld)
             if(msq .gt. maxmsq) cycle  !!ignroe heavy fields
             if(msq .lt. -kbyasq*0.99)then
                write(*,*) "tachyonic field", fld
                write(*,*) "m^2 = ", msq
                write(*,*) "ln a = ", lna
                call coop_return_error("coop_lattice_initial_power_initialize", "not sure how to set initial conditions for tachyonic field(s)", "stop")
             endif
             omega =  sqrt(kbyasq+msq)             
             if(this%lnk(ik) .gt. coop_infbg%lnHend)then
                this%f_cov(fld, fld, ik) = kbya**3/coop_2pi**2/omega
                this%fdbyf(fld, ik) = cmplx(0.d0, -omega)
                cycle
             endif
             call args%init( r = (/ exp(this%lnk(ik)*2.d0) /), i = (/ fld  /) )
             call pert_r%set_arguments(args)
             call pert_i%set_arguments(args)
             ini_r(1) = kbya**1.5/coop_2pi/sqrt(omega)
             ini_r(2) = 0.d0
             ini_i(1) = 0.d0
             ini_i(2) = -omega*ini_r(1)
             call pert_r%set_initial_conditions(xini = lna, yini = ini_r(1:2))
             call pert_i%set_initial_conditions(xini = lna, yini = ini_i(1:2))
             call pert_r%evolve(coop_lattice_perturb_diag_eqs, 0.d0)
             call pert_i%evolve(coop_lattice_perturb_diag_eqs, 0.d0)
             this%f_cov(fld, fld, ik)  = pert_r%y(1) * pert_r%y(1) + pert_i%y(1) * pert_i%y(1)
             this%fdbyf(fld, ik) = cmplx(pert_r%y(2), pert_i%y(2))/cmplx(pert_r%y(1), pert_i%y(1))
          enddo
       enddo
    else
       call pert_r%init(n=2*coop_infbg%nflds, method=COOP_ODE_DVERK, tol = 1.d-8)
       call pert_i%init(n=2*coop_infbg%nflds, method=COOP_ODE_DVERK, tol = 1.d-8)
       do ik = 1, this%nk
          i = 1
          do while(coop_infbg%lna(i)+coop_infbg%lnH(i) .lt. this%lnk(ik) - shift .and. i.lt. coop_infbg%nsteps)
             i = i+1
          enddo
          i_start = i - 1
          lna = coop_infbg%lna(i_start)
          kbya = exp(this%lnk(ik)-lna)
          kbyasq = kbya ** 2
          do fld = 1, coop_infbg%nflds
             msq = coop_infbg%msq(lna, fld, fld)
             if(msq .lt. -kbyasq*0.99)then
                write(*,*) "tachyonic field", fld
                write(*,*) "m^2 = ", msq
                write(*,*) "ln a = ", lna
                call coop_return_error("coop_lattice_initial_power_initialize", "not sure how to set initial conditions for tachyonic field(s)", "stop")
             endif
             omega =  sqrt(kbyasq+msq)             !!for initial condition I only consider diagonal part
             ini_r(fld) = kbya**1.5/coop_2pi/sqrt(omega)
             ini_r(coop_infbg%nflds+fld) = 0.d0
             ini_i(fld) = 0.d0
             ini_i(coop_infbg%nflds+fld) = -omega*ini_r(fld)
          enddo
          call args%init( r = (/ exp(this%lnk(ik)*2.d0) /))
          call pert_r%set_arguments(args)
          call pert_i%set_arguments(args)
          call pert_r%set_initial_conditions(xini = lna, yini = ini_r)
          call pert_i%set_initial_conditions(xini = lna, yini = ini_i)
          call pert_r%evolve(coop_lattice_perturb_eqs, 0.d0)
          call pert_i%evolve(coop_lattice_perturb_eqs, 0.d0)
          do i=1, coop_infbg%nflds
             if(cmplx(pert_r%y(i), pert_i%y(i))  .ne. cmplx(0.d0, 0.d0))then
                this%fdbyf(i, ik) = cmplx(pert_r%y(coop_infbg%nflds+i), pert_i%y(coop_infbg%nflds+i))/cmplx(pert_r%y(i), pert_i%y(i))
             else
                this%fdbyf(i, ik) = cmplx(0.d0, 0.d0)
             endif
             do j=1, coop_infbg%nflds             
                this%f_cov(i, j, ik)  = pert_r%y(i) * pert_r%y(j) + pert_i%y(i) * pert_i%y(j)
             enddo
          enddo
       enddo
    endif
  end subroutine coop_lattice_initial_power_initialize
  

  !!y(1:nflds): \delta\phi
  !!y(nflds+1:2*nflds): \dot{\delta\phi}
  subroutine coop_lattice_perturb_eqs(n, lna, y, yp, args)
    COOP_INT::n, i
    type(coop_arguments)::args
#define KSQ  args%r(1)    
    COOP_REAL::y(n), yp(n), lna, H, msqs(n/2, n/2), kbya2
    H = coop_infbg%Hubble(lna)
    yp(1:n/2) = y(n/2+1:n)/H
    msqs = coop_infbg%msqs(lna)
    kbya2 = KSQ*exp(-2.d0*lna)
    do i=1, n/2
       msqs(i, i) = msqs(i, i) + kbya2 
    enddo
    yp(n/2+1:n) = -3.d0*y(n/2+1:n) - matmul(msqs, y(1:n/2))/H    
#undef KSQ    
  end subroutine coop_lattice_perturb_eqs


  !!args%i(1): i
  !!y(1): \delta\phi_i
  !!y(2): \dot{\delta\phi_i}
  subroutine coop_lattice_perturb_diag_eqs(n, lna, y, yp, args)
    COOP_INT::n, i
    type(coop_arguments)::args
#define KSQ  args%r(1)
#define IND args%i(1)    
    COOP_REAL::y(2), yp(2), lna, H, msq
    H = coop_infbg%Hubble(lna)
    yp(1) = y(2)/H
    msq = coop_infbg%msq(lna, IND, IND)
    yp(2) = -3.d0*y(2) - (KSQ*exp(-2.d0*lna) + msq)*y(1)/H
#undef KSQ
#undef IND    
  end subroutine coop_lattice_perturb_diag_eqs
  
  
end module coop_lattice_fields_mod



