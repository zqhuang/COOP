module coop_ellipse_collapse_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"


  type coop_ellipse_collapse_params
     COOP_REAL::lambda1, lambda2, lambda3 !!lambda's
     COOP_REAL::zinit = 50.d0 !!initial redshift
     COOP_REAL::t = 0.d0
     COOP_INT::num_ode_vars = 6
     type(coop_2dfunction)::bprime
   contains
     procedure::free => coop_ellipse_collapse_params_free
     procedure::init => coop_ellipse_collapse_params_init
  end type coop_ellipse_collapse_params

#define  X_1  y(1)
#define  X_2  y(2)
#define  X_3  y(3)
#define  DX1DT  y(4)
#define  DX2DT  y(5)
#define  DX3DT  y(6)

contains

  subroutine coop_ellipse_collapse_params_free(this)
    class(coop_ellipse_collapse_params)::this
    call this%bprime%free()
  end subroutine coop_ellipse_collapse_params_free


  subroutine coop_ellipse_collapse_params_init(this)
    class(coop_ellipse_collapse_params)::this
    call this%free()
    call this%bprime%init_symmetric(f = coop_ellipse_collapse_bprime_reduced, nx = 200, xmin = 1.d-3, xmax = 1.d3, xlog = .true., name = "BPRIME")
  end subroutine coop_ellipse_collapse_params_init


  
  !!evolve y from t to t_end
  !!vector y = ( x_1, x_2,  x_3,  d x_1/dt, d x_2/dt, d x_3/dt )
  !!the object params contains all the parameters for the model
  subroutine coop_ellipse_collapse_evolve(params, t, y, t_end)
    type(coop_ellipse_collapse_params) params
    COOP_REAL::t, y(params%num_ode_vars), t_end
    COOP_REAL,parameter::tol = 1.d-7
    COOP_INT::ind, nw
    COOP_REAL::c(24), w(params%num_ode_vars, 9)
    ind = 1
    call coop_dverk_with_ellipse_collapse_params(params%num_ode_vars, coop_ellipse_collapse_odes, params, t, y, t_end, tol, ind, c, params%num_ode_vars, w)
  end subroutine coop_ellipse_collapse_evolve



  subroutine coop_ellipse_collapse_odes(n, t, y, dydt, params)
    COOP_INT::n
    COOP_REAL::t, y(n), dydt(n)
    type(coop_ellipse_collapse_params)::params
  end subroutine coop_ellipse_collapse_odes


!!=========================================================================
!!utilities; no need to understand or change anything below
  subroutine coop_dverk_with_ellipse_collapse_params(n, fcn, params, x, y, xend, tol, ind, c, nw, w)
    type(coop_ellipse_collapse_params) params
#define DVERK_ARGUMENTS ,params
#include "dverk.h"    
#undef DVERK_ARGUMENTS
  end subroutine coop_dverk_with_ellipse_collapse_params

!!$  subroutine coop_ellipse_collapse_compute_bprime(a, bprime)
!!$    COOP_REAL::a(3), bprime(3)
!!$    type(coop_arguments)::args
!!$    call args%init(r = (/ (a(1)/a(2))**2,  (a(1)/a(3))**2 /) )
!!$    call calc_int(bprime(1))
!!$    call args%init(r = (/ (a(2)/a(1))**2,  (a(2)/a(3))**2 /) )
!!$    call calc_int(bprime(2))
!!$    call args%init(r = (/ (a(3)/a(1))**2,  (a(3)/a(2))**2 /) )


  function coop_ellipse_collapse_bprime_reduced(lambda1, lambda2) result(bprime)
    COOP_REAL::lambda1, lambda2, eps, bprime
    type(coop_arguments)::args
    call args%init( r = (/ lambda1, lambda2 /) )
    if(abs(lambda1-1.d0) .lt. 0.1d0 .and. abs(lambda2-1.d0) .lt. 0.1d0)then
       bprime = coop_integrate(coop_ellipse_collapse_bprime_int, 0.d0, 1.d0, args, 1.d-8) 
    else
       eps = min(lambda1, lambda2, 0.01d0)/5.d0
       bprime = coop_integrate(coop_ellipse_collapse_b_int, eps, 1.d0, args, 1.d-8) -2.d0/3.d0 + eps*sqrt(eps/args%r(1)/args%r(2)) * ( 2.d0/3.d0+ eps*( 0.4d0*(1.d0-0.5d0/args%r(1)-0.5d0/args%r(2)) + eps*(2.d0/7.d0)*(1.d0-1.d0/args%r(1)-1.d0/args%r(2)+0.25d0/args%r(1)/args%r(2)+0.375d0/args%r(1)**2+0.375d0/args%r(2)**2) ) )
    endif
    call args%free()
  end function coop_ellipse_collapse_bprime_reduced


  function coop_ellipse_collapse_bprime_int(x, args) result(f)
    type(coop_arguments)::args
    COOP_REAL::x, f
    f = sqrt(x)*(1.d0/sqrt((args%r(1)*(1.d0-x)+x)*(args%r(2)*(1.d0-x)+x))-1.d0)
  end function coop_ellipse_collapse_bprime_int

  function coop_ellipse_collapse_b_int(x, args) result(f)
    type(coop_arguments)::args
    COOP_REAL::x, f
    f = sqrt(x / ((args%r(1)*(1.d0-x)+x)*(args%r(2)*(1.d0-x)+x)))
  end function coop_ellipse_collapse_b_int


end module coop_ellipse_collapse_mod
