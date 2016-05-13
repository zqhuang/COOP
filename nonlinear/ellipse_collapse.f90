module coop_ellipse_collapse_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

!!by default everything is a function of scale factor a;  a = 1 today.
!!physical time unit = 1/H_0 

  type coop_ellipse_collapse_params
     COOP_REAL,dimension(3)::lambda = (/ 0.d0, 0.d0, 0.d0 /) !!lambda's
     COOP_REAL::zinit = 100.d0 !!initial redshift
     COOP_REAL::Omega_m = 0.3d0  !!fractional matter density
     COOP_REAL::w = -1.d0   !!dark energy EOS
     COOP_REAL::Omega_r, Omega_de
     COOP_REAL::Omega_k = 0.d0
     COOP_REAL::h = 0.7d0
     COOP_REAL::T_CMB = 2.726  
     logical::is_spherical = .false.
     COOP_INT::num_ode_vars = 6
     COOP_REAL,dimension(3)::collapse_a_ratio = (/ 0.18d0, 0.18d0, 0.18d0 /)
     type(coop_2dfunction)::bprime  
     type(coop_function)::Dbya, Dbyadot
   contains
     procedure::free => coop_ellipse_collapse_params_free  !!subroutine; no argument; release the memory allocated for this object
     procedure::init => coop_ellipse_collapse_params_init  !!subroutine; this%init(Omega_m, w, Omega_k, h, F_pk, e_nu, p_nu);  initialize the object with these parameters
     procedure::get_bprime => coop_ellipse_collapse_params_get_bprime !!subroutine; this%get_bprime(a, bprime)  where a(1:3) is the input array a_1, a_2, a_3, bprime(1:3) is the output b'_1, b'_2, b'_3; 
     procedure::Growth_D => coop_ellipse_collapse_params_Growth_D !! this%Growth_D(a) is a function that returns the growth factor D; input a is the scale factor
     procedure::Growth_H_D => coop_ellipse_collapse_params_Growth_H_D  !!this%Growth_H_D(a) is a function that returns d ln D/d t, where a is the scale factor
     procedure::dadt => coop_ellipse_collapse_params_aH  !!this%dadt(a) is a function that  returns da/dt
     procedure::ddotabya => coop_ellipse_collapse_params_ddotabya !!this%ddotabya is a function that returns ( \ddot a / a )
     procedure::set_initial_conditions => coop_ellipse_collapse_params_set_initial_conditions  !!this%set_initial_conditions(y) is a subroutine set the initial conditions for y = (a_1, a_2, a_3, d a_1/dt, d a_2/dt, d a_3/dt) 
     procedure::evolve=> coop_ellipse_collapse_params_evolve  !!this%evolve(a, y, a_end) is a subroutine that evolves the current  y = (a_1, a_2, a_3, d a_1/dt, d a_2/dt, d a_3/dt) from current scale factor a to the scale factor a_end; the input is y at a; after the call a becomes a_end.
     procedure::set_growth_initial_conditions => coop_ellipse_collapse_params_set_growth_initial_conditions !!set initial conditions for the growth function solver
  end type coop_ellipse_collapse_params


contains

  subroutine coop_ellipse_collapse_params_set_initial_conditions(this, y)
    !!set inital vector y = ( a_1, a_2,  a_3,  d a_1/dt, d a_2/dt, d a_3/dt ) at redshift this%zinit
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::y(:), D_ini, a_ini, dadt_ini
    a_ini = 1.d0/(1.d0+this%zinit)
    D_ini = this%Growth_D(a_ini)
    dadt_ini = this%dadt(a_ini)
    y(1:3) = a_ini *(1.d0-this%lambda*D_ini)
    y(4:6) = y(1:3)*dadt_ini/a_ini - a_ini*this%lambda*D_ini*this%Growth_H_D(a_ini)
  end subroutine coop_ellipse_collapse_params_set_initial_conditions

  subroutine coop_ellipse_collapse_odes(n, a, y, dyda, params)
  !!the ODE that evolves  y = ( a_1, a_2,  a_3,  d a_1/dt, d a_2/dt, d a_3/dt )
  !!other inputs: 
  !!n = 6 is the dimension of y
  !!a is the scale factor
  !!params is the object contain all the parameters and methods
  !! return dyda = d y/d a
    COOP_INT::n
    COOP_REAL::a, y(n), dyda(n), dadt, bprime(3), growthD, delta, dark_Energy_term, rhomby3, radiation_term
    type(coop_ellipse_collapse_params)::params
    dadt = params%dadt(a)
    delta = a**3/(y(1)*y(2)*y(3))-1.d0
    radiation_term = - params%omega_r/a**4*2.d0
    dark_Energy_term =  - params%omega_de*a**(-3.d0*(1.d0+params%w))*(1.d0+3.d0*params%w)  !!dark energy contribution; ignore dark energy perturbations in wCDM
    rhomby3 = params%Omega_m/a**3 !!I am working in unit of H_0^2/(8\pi G)
    if(params%is_spherical)then
       if(y(1)/a/params%collapse_a_ratio(1).lt. 0.9999d0 .and. y(4) .lt. 0.d0)then  !!collapsed; freeze it
          dyda(1) = 0.d0
          dyda(4) = 0.d0
       else  !!still collapsing
          dyda(1) = y(4) / dadt    !!d a_1/da = d a_1/dt /(da/dt)
          dyda(4) = y(1)/dadt/2.d0 * ( &   !! d( da_1/dt)/da =(d^2 a_1/dt^2)/(da/dt)
               dark_Energy_term  &
               + radiation_term &
               -  rhomby3*(1.d0 + delta) &  !!matter contribution
               )
       endif
       dyda(2:3) = dyda(1)
       dyda(5:6) = dyda(4)
    else
       call params%get_bprime(y(1:3), bprime)
       growthD = params%growth_D(a)
       if(y(1)/a/params%collapse_a_ratio(1).lt. 0.9999d0 .and. y(4) .lt. 0.d0)then  !!collapsed; freeze it
          dyda(1) = 0.d0
          dyda(4) = 0.d0
       else  !!still collapsing
          dyda(1) = y(4) / dadt    !!d a_1/da = d a_1/dt /(da/dt)
          dyda(4) = y(1)/dadt/2.d0 * ( &   !! d( da_1/dt)/da =(d^2 a_1/dt^2)/(da/dt)
               dark_Energy_term  &
               + radiation_term &
               -  rhomby3*(1.d0 + delta *(1.d0+bprime(1)*1.5d0) + (3.d0*params%lambda(1)-sum(params%lambda))*growthD ) &  !!matter contribution
               )
       endif
       if(y(2)/a /params%collapse_a_ratio(2).lt. 0.9999d0 .and. y(5) .lt. 0.d0 )then !!collapsed; freeze it
          dyda(2) = 0.d0
          dyda(5) = 0.d0
       else !!still collapsing
          dyda(2) = y(5) / dadt    !!d a_1/da = d a_1/dt /(da/dt)
          dyda(5) = y(2)/dadt/2.d0 * ( &   !! d( da_1/dt)/da =(d^2 a_1/dt^2)/(da/dt)
               dark_Energy_term  &
               + radiation_term &  
               - rhomby3 *(1.d0 + delta *(1.d0+bprime(2)*1.5d0) + (3.d0*params%lambda(2)-sum(params%lambda))*growthD ) &  !!matter contribution
               )
       endif
       if(y(3)/a /params%collapse_a_ratio(3).lt. 0.9999d0.and. y(6) .lt. 0.d0)then !!collapsed; freeze it
          dyda(3) = 0.d0
          dyda(6) = 0.d0
       else !!still collapsing
          dyda(3) = y(6) / dadt    !!d a_1/da = d a_1/dt /(da/dt)
          dyda(6) = y(3)/dadt/2.d0 * ( &   !! d( da_1/dt)/da =(d^2 a_1/dt^2)/(da/dt)
               dark_Energy_term  &
               + radiation_term &  
               -  rhomby3*(1.d0 + delta *(1.d0+bprime(3)*1.5d0) + (3.d0*params%lambda(3)-sum(params%lambda))*growthD ) &  !!matter contribution
               )       
       endif
    endif
  end subroutine coop_ellipse_collapse_odes


!!y = ( D/a,  d (D/a)/dt )
!!set initial conditions at this%zinit
  subroutine coop_ellipse_collapse_params_set_growth_initial_conditions(this, y)
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a, dadt, addbya, H
    COOP_REAL::y(2)
    a = 1.d0/(1.d0+this%zinit)
    dadt = this%dadt(a)
    H = dadt/a
    addbya = this%ddotabya(a)
    y(1) = 1.d0
    y(2) =  - (2.d0*H**2+addbya-1.5d0*this%Omega_m/a**3)*y(1)/4.d0/H
  end subroutine coop_ellipse_collapse_params_set_growth_initial_conditions


!!y = ( D/a,  d (D/a)/dt )
!!the equation in the standard wCDM case is
!!\ddot (aD) + 4 H (aD) + (2H^2 - \ddot a/ a - 3/2 Omega_m0 / a^3) (aD) = 0
!!n is the dimension of y (here = 2)
!!a is the scale factor
!!dyda returns dy/da
!!params is the object that contains all the parameters
  subroutine coop_ellipse_collapse_growth_ode(n, a, y, dyda, params)
    type(coop_ellipse_collapse_params)::params
    COOP_INT::n
    COOP_REAL::a, dadt, addbya, H
    COOP_REAL::y(2), dyda(2)
    dadt = params%dadt(a)
    H = dadt/a
    addbya = params%ddotabya(a)
    dyda(1) = y(2)/dadt  !! d(D/a)/ da  = d(D/a)/dt / (da/dt)
    dyda(2) = (-4.d0*H*y(2) - (2.d0*H**2+addbya-1.5d0*params%Omega_m/a**3)*y(1))/dadt
  end subroutine coop_ellipse_collapse_growth_ode



!!=====================You don't need to read anything below ==================

!!==================== methods of class coop_ellipse_collapse_params ============
  subroutine coop_ellipse_collapse_params_evolve(this, a, y, a_end)
  !!evolve y from a to a_end
  !!vector y = ( a_1, a_2,  a_3,  d a_1/dt, d a_2/dt, d a_3/dt )
  !!the object this contains all the parameters for the model
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a, y(this%num_ode_vars), a_end
    COOP_REAL,parameter::tol = 1.d-8
    COOP_INT::ind
    COOP_REAL::c(24), w(this%num_ode_vars, 9)
    select type(this)
    type is(coop_ellipse_collapse_params)
       ind = 1
       call coop_dverk_with_ellipse_collapse_params(this%num_ode_vars, coop_ellipse_collapse_odes, this, a, y, a_end, tol, ind, c, this%num_ode_vars, w)
    class default
       stop "Evolve: Extended class this has not been implemented"
    end select
  end subroutine coop_ellipse_collapse_params_evolve


  subroutine coop_ellipse_collapse_params_init(this, Omega_m, w, h, Omega_k, lambda, F_pk, e_nu, p_nu)
    class(coop_ellipse_collapse_params)::this
    COOP_INT,parameter::na = 256
    COOP_REAL::a(na), Dbya(na), Dbyadot(na), adynamic
    COOP_INT::i
    COOP_REAL,optional::Omega_m, w, Omega_k, h, lambda(3), F_pk, e_nu, p_nu
    logical cosmology_updated
    !!!!for D(a) solver
    COOP_REAL::y(2)  
    COOP_REAL,parameter::tol = 1.d-8
    COOP_INT::ind
    COOP_REAL::c(24), wspace(2, 9)

    cosmology_updated = .false.
    if(present(Omega_m))then
       if(abs(this%Omega_m - Omega_m) .gt. 1.d-5)then
          this%Omega_m = Omega_m
          cosmology_updated = .true.
       endif
    endif
    if(present(Omega_k))then
       if(abs(this%Omega_k - Omega_k) .gt. 1.d-5)then
          this%Omega_k = Omega_k
          cosmology_updated = .true.
       endif
    endif
    if(present(w))then
       if(abs(this%w - w) .gt. 1.d-5)then
          this%w = w
          cosmology_updated = .true.
       endif
    endif
    if(present(h))then
       if(abs(this%h - h) .gt. 1.d-5)then
          this%h = h
          cosmology_updated = .true.
       endif
    endif
    this%Omega_r = 4.187d-5/this%h**2*(this%T_CMB/2.726)**4
    this%Omega_de = 1.d0 - this%Omega_m - this%Omega_r - this%Omega_k
    if(present(lambda))then
       this%lambda = lambda
       if(present(F_pk) .or. present(p_nu) .or. present(e_nu))then
          stop "You pass either lambda or (F_pk, e_nu, p_nu) to init, not both."
       endif
    elseif(present(F_pk) .and. present(p_nu) .and. present(e_nu))then
       if(F_pk .lt. 0.d0 .or. e_nu .lt. 0.d0 .or. abs(p_nu) .gt. e_nu)then
          write(*,*) "Invalid F_pk, e_nu, p_nu:"
          write(*,*) F_pk, e_nu, p_nu
          write(*,*) "check conditions: F_pk >=0, e_nu >= 0, and |p_nu| <= e_nu"
          stop
       endif
       this%lambda(3) = (F_pk/3.d0)*(1.d0 + 3.d0*e_nu + p_nu)
       this%lambda(2) = (F_pk/3.d0)*(1.d0 - 2.d0*p_nu)
       this%lambda(1) = (F_pk/3.d0)*(1.d0 - 3.d0*e_nu + p_nu)
    endif
    this%is_spherical = abs(this%lambda(1) - this%lambda(2)) .lt. 1.d-8 .and. abs(this%lambda(1) - this%lambda(3)) .lt. 1.d-8 .and. abs(this%collapse_a_ratio(1) - this%collapse_a_ratio(2)) .lt. 1.d-3 .and. abs(this%collapse_a_ratio(1)-this%collapse_a_ratio(3)) .lt. 1.d-3
    !!set b' function
    if(.not. this%bprime%initialized)call this%bprime%init_symmetric(f = coop_ellipse_collapse_bprime_reduced, nx = 301, xmin = 1.d-5, xmax = 1.d5, xlog = .true., name = "BPRIME")

    !!set  D(a)/a
    select type(this)
    type is(coop_ellipse_collapse_params)
       
       if( cosmology_updated .or. .not.this%Dbya%initialized)then
          call this%Dbya%free()
          call this%Dbyadot%free()
          call coop_set_uniform(na, a, 1.d0/(1.d0+ this%zinit), 1.d0)
          call this%set_Growth_initial_conditions(y)
          Dbya(1) = y(1)
          Dbyadot(1) = y(2)
          adynamic = a(1)
          ind = 1
          do i  = 2, na
             call coop_dverk_with_ellipse_collapse_params(2, coop_ellipse_collapse_growth_ode, this, adynamic, y, a(i), tol, ind, c, 2, wspace)
             Dbya(i) = y(1)
             Dbyadot(i) = y(2)
          enddo
          dbyadot = dbyadot/dbya(na)
          dbya = dbya/dbya(na)
          call this%Dbya%init(n = na, xmin = 1.d0/(1.d0+this%zinit), xmax = 1.d0, f = Dbya, method=COOP_INTERPOLATE_SPLINE, check_boundary = .false., name="D(a)/a")
          call this%Dbyadot%init(n = na, xmin = 1.d0/(1.d0+this%zinit), xmax = 1.d0, f = Dbyadot, method=COOP_INTERPOLATE_SPLINE, check_boundary = .false., name="dot(D(a)/a)")
       endif
    class default
       stop "extended class for ode is not yet supported"
    end select
  end subroutine coop_ellipse_collapse_params_init


  subroutine coop_ellipse_collapse_params_get_bprime(this, a, bprime)
    class(coop_ellipse_collapse_params)::this
    COOP_REAL,intent(IN)::a(3)
    COOP_REAL,intent(OUT)::bprime(3)
    bprime(1) = this%bprime%eval( (a(1)/a(2))**2,  (a(1)/a(3))**2 )
    bprime(2) = this%bprime%eval( (a(2)/a(1))**2,  (a(2)/a(3))**2 )
    bprime(3) = this%bprime%eval( (a(3)/a(2))**2,  (a(3)/a(1))**2 )
  end subroutine coop_ellipse_collapse_params_get_bprime

  !!D(a)
  function coop_ellipse_collapse_params_Growth_D(this, a) result(D)
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a, D
    D = this%Dbya%eval(a)*a
  end function coop_ellipse_collapse_params_Growth_D

  !!d ln D/dt 
  function coop_ellipse_collapse_params_Growth_H_D(this, a) result(H_D)
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a, H_D
    H_D = this%Dbyadot%eval(a)/this%Dbya%eval(a) + this%dadt(a)/a
  end function coop_ellipse_collapse_params_Growth_H_D
  
  !! H a / (H_0 a_0)
  function coop_ellipse_collapse_params_aH(this, a) result(aH)
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a, aH
    aH = sqrt(this%Omega_k + (this%Omega_m+this%Omega_r/a)/a + this%omega_de*a**(-1.d0-3.d0*this%w))
  end function coop_ellipse_collapse_params_aH

!!return H_0^2  \ddot a / a
  function coop_ellipse_collapse_params_ddotabya(this, a) result(ddotabya)
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a, ddotabya
    ddotabya = -( (this%Omega_r/a*2.d0 + this%Omega_m)/a**3 + this%Omega_de*a**(-3.d0*(1.d0+this%w))*(1.d0+3.d0*this%w) )/2.d0
  end function coop_ellipse_collapse_params_ddotabya

  subroutine coop_ellipse_collapse_params_free(this)
    class(coop_ellipse_collapse_params)::this
    call this%bprime%free()
    call this%Dbya%free()
    call this%Dbyadot%free()
  end subroutine coop_ellipse_collapse_params_free


!!=========================================================================
!!utilities; no need to understand or change anything below
  subroutine coop_dverk_with_ellipse_collapse_params(n, fcn, params, x, y, xend, tol, ind, c, nw, w)
    type(coop_ellipse_collapse_params) params
#define DVERK_ARGUMENTS ,params
#include "dverk.h"    
#undef DVERK_ARGUMENTS
  end subroutine coop_dverk_with_ellipse_collapse_params


  function coop_ellipse_collapse_bprime_reduced(lambda1, lambda2) result(bprime)
    COOP_REAL::lambda1, lambda2, eps, bprime, xcut
    type(coop_arguments)::args
    if(abs(lambda1-1.d0) .lt. 1.d-3 .and. abs(lambda2-1.d0) .lt. 1.d-3)then !!both close to 0
       bprime = (2.d0/15.d0)*(2.d0-lambda1 - lambda2) + (3.d0/28.d0-0.05d0)*((1.d0-lambda1)**2+(1.d0-lambda2)**2+(2.d0/3.d0)*(1.d0-lambda1)*(1.d0-lambda2))
       return
    endif
    call args%init( r = (/ lambda1, lambda2 /) )
    if(max(lambda1, lambda2) .gt. 1.d2 .or. lambda1*lambda2 .gt. 1.d0)then
       xcut = min(1.d3/max(lambda1, lambda2)**0.25, 1.d2/(lambda1*lambda2)**(1.d0/6.d0))
       bprime = (-2.d0/3.d0) &
         + 2.d0*coop_integrate(coop_ellipse_collapse_b_int_form3, 0.d0, xcut, args, 1.d-6)
       return
    endif
    xcut = 100.d0
    bprime = (-2.d0/3.d0)+ coop_integrate(coop_ellipse_collapse_b_int_form1, 0.d0, 1.d0/(1.d0+xcut), args, 1.d-6)  &
         + coop_integrate(coop_ellipse_collapse_b_int_form2, 0.d0, xcut, args, 1.d-6)
    call args%free()
  end function coop_ellipse_collapse_bprime_reduced

  function coop_ellipse_collapse_bprime_accurate(lambda1, lambda2) result(bprime)
    COOP_REAL::lambda1, lambda2, eps, bprime, xcut
    type(coop_arguments)::args
    if(abs(lambda1-1.d0) .lt. 1.d-4 .and. abs(lambda2-1.d0) .lt. 1.d-4)then !!both close to 0
       bprime = (2.d0/15.d0)*(2.d0-lambda1 - lambda2) + (3.d0/28.d0-0.05d0)*((1.d0-lambda1)**2+(1.d0-lambda2)**2+(2.d0/3.d0)*(1.d0-lambda1)*(1.d0-lambda2))
       return
    endif
    call args%init( r = (/ lambda1, lambda2 /) )
    if(max(lambda1, lambda2) .gt. 1.d2 .or. lambda1*lambda2 .gt. 1.d0)then
       xcut = min(1.d4/max(lambda1, lambda2)**0.25, 1.d3/(lambda1*lambda2)**(1.d0/6.d0))
       bprime = (-2.d0/3.d0) &
         + 2.d0*coop_integrate(coop_ellipse_collapse_b_int_form3, 0.d0, xcut, args, 1.d-8)
       return
    endif
    xcut = 100.d0
    bprime = (-2.d0/3.d0)+ coop_integrate(coop_ellipse_collapse_b_int_form1, 0.d0, 1.d0/(1.d0+xcut), args, 1.d-8)  &
         + coop_integrate(coop_ellipse_collapse_b_int_form2, 0.d0, xcut, args, 1.d-8)
    call args%free()
  end function coop_ellipse_collapse_bprime_accurate


  function coop_ellipse_collapse_b_int_form1(x, args) result(f)
    type(coop_arguments)::args
    COOP_REAL::x, f
    f = sqrt(x / ((args%r(1)*(1.d0-x)+x)*(args%r(2)*(1.d0-x)+x)))
  end function coop_ellipse_collapse_b_int_form1


  function coop_ellipse_collapse_b_int_form2(x, args) result(f)
    type(coop_arguments)::args
    COOP_REAL::x, f
    f = 1.d0/(1.d0+x)/sqrt((1.d0+x)*(1.d0+args%r(1)*x)*(1.d0+args%r(2)*x))
  end function coop_ellipse_collapse_b_int_form2

  function coop_ellipse_collapse_b_int_form3(x, args) result(f)
    type(coop_arguments)::args
    COOP_REAL::x, f
    f = x/(1.d0+x**2)/sqrt((1.d0+x**2)*(1.d0+args%r(1)*x**2)*(1.d0+args%r(2)*x**2))
  end function coop_ellipse_collapse_b_int_form3




end module coop_ellipse_collapse_mod
