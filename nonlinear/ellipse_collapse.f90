module coop_ellipse_collapse_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

!!by default everything is a function of scale factor a;  a = 1 today.
!!physical time unit = 1/H_0 

  type coop_ellipse_collapse_params
     COOP_REAL,dimension(3)::lambda = (/ 0.d0, 0.d0, 0.d0 /) !!lambda's
     COOP_REAL::zinit = 50.d0 !!initial redshift
     COOP_REAL::Omega_m = 0.3d0  !!fractional matter density
     COOP_REAL::w = -1.d0   !!dark energy EOS
     COOP_REAL::t = 0.d0
     COOP_INT::num_ode_vars = 6
     COOP_REAL,dimension(3)::collapse_a_ratio = (/ 0.18d0, 0.18d0, 0.18d0 /)
     type(coop_2dfunction)::bprime  !!tabulated b' function
     type(coop_function)::Dbya   !!growth factor D(a)/a = delta (a) / a delta(today); its asymptotic limit -> constant when a<<1;
   contains
     procedure::free => coop_ellipse_collapse_params_free
     procedure::init => coop_ellipse_collapse_params_init
     procedure::get_bprime => coop_ellipse_collapse_params_get_bprime
     procedure::Growth_D => coop_ellipse_collapse_params_Growth_D !! D(a)
     procedure::Growth_H_D => coop_ellipse_collapse_params_Growth_H_D  !!d ln D/d t as a function of a
     procedure::delta_m => coop_ellipse_collapse_params_delta_m
     procedure::dadt => coop_ellipse_collapse_params_aH
     procedure::set_initial_conditions => coop_ellipse_collapse_params_set_initial_conditions
     procedure::evolve=> coop_ellipse_collapse_params_evolve
  end type coop_ellipse_collapse_params


contains

  subroutine coop_ellipse_collapse_params_evolve(this, a, y, a_end)
  !!evolve y from a to a_end
  !!vector y = ( a_1, a_2,  a_3,  d a_1/dt, d a_2/dt, d a_3/dt )
  !!the object this contains all the parameters for the model
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a, y(this%num_ode_vars), a_end
    COOP_REAL,parameter::tol = 1.d-7
    COOP_INT::ind, nw
    COOP_REAL::c(24), w(this%num_ode_vars, 9)
    select type(this)
    type is(coop_ellipse_collapse_params)
       ind = 1
       call coop_dverk_with_ellipse_collapse_params(this%num_ode_vars, coop_ellipse_collapse_odes, this, a, y, a_end, tol, ind, c, this%num_ode_vars, w)
    class default
       stop "Evolve: Extended class this has not been implemented"
    end select
  end subroutine coop_ellipse_collapse_params_evolve

  subroutine coop_ellipse_collapse_params_set_initial_conditions(this, y)
  !!set inital vector y = ( a_1, a_2,  a_3,  d a_1/dt, d a_2/dt, d a_3/dt )
  !!this is the object contain all the parameters and methods
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::y(:), D_ini, a_ini, dadt_ini
    a_ini = 1.d0/(1.d0+this%zinit)
    D_ini = this%Growth_D(a_ini)
    dadt_ini = this%dadt(a_ini)
    y(1:3) = a_ini *(1.d0-this%lambda*D_ini)
    y(4:6) = y(1:3)*dadt_ini/a_ini - a_ini*this%lambda*this%Growth_H_D(a_ini)
  end subroutine coop_ellipse_collapse_params_set_initial_conditions

  subroutine coop_ellipse_collapse_odes(n, a, y, dyda, params)
  !!the main ODE
  !!input 
  !!vector y = ( a_1, a_2,  a_3,  d a_1/dt, d a_2/dt, d a_3/dt )
  !!n = 6 is the dimension of y
  !!a is the scale factor
  !!params is the object contain all the parameters and methods
  !! return dyda = d y/d a
    COOP_INT::n
    COOP_REAL::a, y(n), dyda(n), dadt, bprime(3), growthD, delta, dark_Energy_term, rhom
    type(coop_ellipse_collapse_params)::params
    dadt = params%dadt(a)
    call params%get_bprime(y(1:3), bprime)
    growthD = params%growth_D(a)
    delta = a**3/(y(1)*y(2)*y(3))-1.d0
    dark_Energy_term =  -(1.d0-params%Omega_m)*a**(-3.d0*(1.d0+params%w))*(1.d0+3.d0*params%w)  !!dark energy contribution; ignore dark energy perturbations in wCDM
    rhom = params%Omega_m/a**3
    if(y(1)/a .le. params%collapse_a_ratio(1))then  !!collapsed; freeze it
       dyda(1) = 0.d0
       dyda(4) = 0.d0
    else  !!still collapsing
       dyda(1) = y(4) / dadt    !!d a_1/da = d a_1/dt /(da/dt)
       dyda(4) = y(1)/dadt/2.d0 * ( &   !! d( da_1/dt)/da =(d^2 a_1/dt^2)/(da/dt)
            dark_Energy_term  &  
            -  rhom*(1.d0 + delta *(1.d0+bprime(1)*1.5d0) + params%lambda(1)*growthD ) &  !!matter contribution
            )
    endif
    if(y(2)/a .le. params%collapse_a_ratio(2))then !!collapsed; freeze it
       dyda(2) = 0.d0
       dyda(5) = 0.d0
    else !!still collapsing
       dyda(2) = y(5) / dadt    !!d a_1/da = d a_1/dt /(da/dt)
       dyda(5) = y(2)/dadt/2.d0 * ( &   !! d( da_1/dt)/da =(d^2 a_1/dt^2)/(da/dt)
            dark_Energy_term  &  
            - rhom *(1.d0 + delta *(1.d0+bprime(2)*1.5d0) + params%lambda(2)*growthD ) &  !!matter contribution
            )
    endif
    if(y(3)/a .le. params%collapse_a_ratio(3))then !!collapsed; freeze it
       dyda(3) = 0.d0
       dyda(6) = 0.d0
    else !!still collapsing
       dyda(3) = y(6) / dadt    !!d a_1/da = d a_1/dt /(da/dt)
       dyda(6) = y(3)/dadt/2.d0 * ( &   !! d( da_1/dt)/da =(d^2 a_1/dt^2)/(da/dt)
            dark_Energy_term  &  
            -  rhom*(1.d0 + delta *(1.d0+bprime(3)*1.5d0) + params%lambda(3)*growthD ) &  !!matter contribution
            )       
    endif
  end subroutine coop_ellipse_collapse_odes
!!==================== methods of class coop_ellipse_collapse_params ============
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
    H_D = (this%Dbya%derivative(a)*a**2/this%Growth_D(a)+1.d0)*this%dadt(a)/a
  end function coop_ellipse_collapse_params_Growth_H_D
  
  !! H a / (H_0 a_0)
  function coop_ellipse_collapse_params_aH(this, a) result(aH)
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a, aH
    aH = sqrt(this%Omega_m/a + (1.d0-this%Omega_m)*a**(-1.d0-3.d0*this%w))
  end function coop_ellipse_collapse_params_aH

  function coop_ellipse_collapse_params_delta_m(this, a) result(delta_m)
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a, delta_m
    delta_m = sum(this%lambda)*this%Growth_D(a)
  end function coop_ellipse_collapse_params_delta_m

  subroutine coop_ellipse_collapse_params_free(this)
    class(coop_ellipse_collapse_params)::this
    call this%bprime%free()
    call this%Dbya%free()
  end subroutine coop_ellipse_collapse_params_free


  subroutine coop_ellipse_collapse_params_init(this, Omega_m, w,  lambda, F_pk, e_nu, p_nu)
    class(coop_ellipse_collapse_params)::this
    COOP_INT,parameter::na = 1024
    COOP_REAL::a(na), Dbya(na)
    COOP_INT::i
    COOP_REAL,optional::Omega_m, w, lambda(3), F_pk, e_nu, p_nu
    logical cosmology_updated
    call this%free()
    cosmology_updated = .false.
    if(present(Omega_m))then
       if(abs(this%Omega_m - Omega_m) .gt. 1.d-5)then
          this%Omega_m = Omega_m
          cosmology_updated = .true.
       endif
    endif
    if(present(w))then
       if(abs(this%w - w) .gt. 1.d-5)then
          this%w = w
          cosmology_updated = .true.
       endif
    endif
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

    !!set b' function
    if(.not. this%bprime%initialized)call this%bprime%init_symmetric(f = coop_ellipse_collapse_bprime_reduced, nx = 301, xmin = 1.d-5, xmax = 1.d5, xlog = .true., name = "BPRIME")

    !!set  D(a)/a
    if( cosmology_updated .or. .not.this%Dbya%initialized)then
       call coop_set_uniform(na, a, 1.d0/(1.d0+ this%zinit), 1.d0)
       !$omp parallel do
       do i  = 1, na
          Dbya(i) = coop_Growth_fitting(Omega_m = this%Omega_m, w = this%w, z = 1.d0/a(i)-1.d0)/a(i)
       enddo
       !$omp end parallel do
       call this%Dbya%init(n = na, xmin = 1.d0/(1.d0+this%zinit), xmax = 1.d0, f = Dbya, method=COOP_INTERPOLATE_SPLINE, check_boundary = .false., name="D(a)/a")
    endif
  end subroutine coop_ellipse_collapse_params_init


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
