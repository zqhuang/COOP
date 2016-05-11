module coop_ellipse_collapse

  type coop_ellipse_collapse_params
     COOP_REAL::lambda1, lambda2, lambda3 !!lambda's
     COOP_REAL::zinit !!initial redshift
     COOP_REAL::t = 0.d0
     COOP_INT::num_ode_vars = 6
  end type coop_ellipse_collapse_params

#define  X_1  y(1)
#define  X_2  y(2)
#define  X_3  y(3)
#define  DX1DT  y(4)
#define  DX2DT  y(5)
#define  DX3DT  y(6)

contains

  
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



  subroutine coop_ellipse_collapse_odes(t, y, dydt)
    
  end subroutine coop_ellipse_collapse_odes


!!=========================================================================
!!utilities; no need to understand or change anything below
  subroutine coop_dverk_with_ellipse_collapse_params(n, fcn, params, x, y, xend, tol, ind, c, nw, w)
    type(coop_ellipse_collapse_params) params
#define DVERK_ARGUMENTS ,params
#include "dverk.h"    
#undef DVERK_ARGUMENTS
  end subroutine coop_dverk_with_ellipse_collapse_params

  subroutine coop_ellipse_collapse_compute_bprime(a, bprime)
    COOP_REAL::a(3), bprime(3)
    COOP_REAL,parameter::xcut = 100.d0
    type(coop_arguments)::args
    call args%init(r = (/ (a(1)/a(2))**2,  (a(1)/a(3))**2 /)
    bprime(1) = coop_integrate(coop_ellipse_collapse_b_int, 0.d0, xcut) + (2.d0/3.d0)/sqrt((1.d0+xcut)**3)*(a(1)**2/a(2)/a(3)) - (2.d0/3.d0)
    call args%init(r = (/ (a(2)/a(1))**2,  (a(2)/a(3))**2 /)
    bprime(2) = coop_integrate(coop_ellipse_collapse_b_int, 0.d0, xcut) + (2.d0/3.d0)/sqrt((1.d0+xcut)**3)*(a(2)**2/a(1)/a(3)) - (2.d0/3.d0)
    call args%init(r = (/ (a(3)/a(1))**2,  (a(3)/a(2))**2 /)
    bprime(2) = coop_integrate(coop_ellipse_collapse_b_int, 0.d0, xcut) + (2.d0/3.d0)/sqrt((1.d0+xcut)**3)*(a(3)**2/a(1)/a(2)) - (2.d0/3.d0)
  end subroutine coop_ellipse_collapse_compute_bprime

  function coop_ellipse_collapse_b_int(x, args) result(f)
    type(coop_arguments)::args
    COOP_REAL::x, f
    f  = 1.d0/(1.d0+x)/sqrt((1.d0+x)*(1.d0+x*args%r(1))*(1.d0+x*args%r(2)))
  end function coop_ellipse_collapse_b_int

end module coop_ellipse_collapse
