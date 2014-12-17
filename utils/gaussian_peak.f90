module coop_gaussian_peak_stat_mod
  use coop_wrapper_typedef
  use coop_special_function_mod
  use coop_integrate_mod
  implicit none
#include  "constants.h"

#define GF_DIM  args%i(1)
#define GF_SIGMA0  args%r(1)
#define GF_SIGMA1 args%r(2)
#define GF_SIGMA2  args%r(3)
#define GF_SIN_BETA args%r(4)
#define GF_COS_BETA args%r(5)
#define GF_CONST args%r(6)    

  

contains

  subroutine coop_gaussian_npeak_set_args(args, dim, sigma0, sigma1, sigma2)
    COOP_INT dim
    COOP_REAL sigma0, sigma1, sigma2, gam
    type(coop_arguments)::args
    gam = sigma0 * sigma2
    if(gam .le. sigma1**2)call coop_return_error("coop_Gaussian_npeak_set_args", "gamma>=1", "stop")
    gam = sigma1**2/gam
    call args%init( i = (/ dim /), r = (/ sigma0, sigma1, sigma2, sqrt(1.d0-gam**2), gam, (sigma2/(sqrt(dble(dim))*sigma1))**dim/coop_2pi**(dim/2.d0+1) /) )
  end subroutine coop_gaussian_npeak_set_args

  function  coop_gaussian_peak_f(v, args) result(f)
    COOP_REAL v, f
    COOP_REAL,parameter::sqrt5h = sqrt(2.5d0)
    type(coop_arguments)::args
    select case(GF_DIM)
    case(1)
       f = abs(v)
    case(2)
       f = v**2-1.d0 + exp(-v**2)
    case(3)
       f = (3.d0 - v**2)*v*(erf(-sqrt5h * v) + erf(-sqrt5h/2.d0 * v))/2.d0 + (1.d0/sqrt5h/coop_sqrtpi)*((31.d0/4.d0 * v**2 + 8.d0/5.d0)*exp(-5.d0/8.d0*v**2) + (v**2/2.d0-8.d0/5.d0)*exp(-2.5d0*v**2))
    case default
       stop "coop_gaussian_peak_f only support dim = 1, 2, 3"
    end select
  end function coop_gaussian_peak_f
    
  !!u, v, dnpk
  function coop_gaussian_npeak_differential(u, v, args) result(dnpk)
    COOP_REAL u, v, dnpk
    type(coop_arguments)::args
    dnpk = GF_CONST/GF_SIN_BETA * coop_gaussian_peak_f(v, args) * exp(- (u**2 + ((v + u*GF_COS_BETA)/GF_SIN_BETA)**2)/2.d0)
  end function coop_gaussian_npeak_differential


  function coop_gaussian_nmax_differential(u, args) result(dnpk)
    COOP_REAL u, v, dnpk
    type(coop_arguments)::args
    dnpk = GF_CONST * coop_gaussian_peak_G(u, args) 
  end function coop_gaussian_nmax_differential

  function coop_gaussian_peak_G(u, args) result(g)
    COOP_REAL u, g, t, q
    type(coop_arguments)::args    
    select case(GF_DIM)
    case(1)
       t = u * GF_COS_BETA/GF_SIN_BETA/coop_sqrt2       
       g = exp(- t**2)*GF_SIN_BETA + sqrt(coop_pi/2.d0)*u * GF_COS_BETA*(1.d0 + erf(t))
    case(2)
       t = u * GF_COS_BETA/GF_SIN_BETA/coop_sqrt2       
       q = sqrt(3.d0-2.d0*GF_COS_BETA**2)
       g = sqrt(coop_pio2)/q*(1.d0+erf(t/q))*exp(-(u*GF_COS_BETA/q)**2) &
            + sqrt(coop_pio2)*(u**2-1.d0)*GF_COS_BETA**2*(1.d0+erf(t)) &
            + u*GF_SIN_BETA*GF_COS_BETA*exp(-t**2)
    case(3)
       t = -u*GF_COS_BETA
       q = 5.d0*GF_SIN_BETA
       if(t + q .lt. 0.d0 .and. t .lt. -3.d0)then
          g = sqrt(coop_2pi)*(3.d0*GF_COS_BETA**2 - t**2)*t
       else
          g = coop_integrate(dgdv, min(t, 0.d0)-q, min(0.d0, t+q), 1.d-7) / GF_SIN_BETA
       endif
    case default
       stop "gaussian_peak_G only support dim = 1, 2, 3"
    end select
    g= g*exp(-u**2/2.d0)
  contains
    function dGdv(v)
      COOP_REAL dGdv, v
      dGdv =  coop_gaussian_peak_f(v, args) * exp(- (((v + u*GF_COS_BETA)/GF_SIN_BETA)**2)/2.d0)
    end function dGdv
  end function coop_gaussian_peak_G

  function coop_gaussian_peak_G_slow(u, args) result(g)
    COOP_REAL t, q, g, u
    type(coop_arguments)::args    
    t = -u*GF_COS_BETA
    q = 6.d0*GF_SIN_BETA
    g = coop_integrate(dgdv, min(t, 0.d0)-q, min(0.d0, t+q), 1.d-8) / GF_SIN_BETA*exp(-u**2/2.d0)
  contains
    function dGdv(v)
      COOP_REAL dGdv, v
      dGdv =  coop_gaussian_peak_f(v, args) * exp(- (((v + u*GF_COS_BETA)/GF_SIN_BETA)**2)/2.d0)
    end function dGdv
  end function coop_gaussian_peak_G_slow


  function coop_gaussian_peak_Gv(u, args) result(g)
    COOP_REAL u, g, t, q
    type(coop_arguments)::args
    
    select case(GF_DIM)
    case(1)
       q = u * GF_COS_BETA       
       t = q/GF_SIN_BETA/coop_sqrt2
       g =  - sqrt(coop_pio2)*(q**2 + GF_SIN_BETA**2)*(1.d0 + Erf(t)) &
            - exp(-t**2)*q*GF_SIN_BETA
    case(2)
       t = u * GF_COS_BETA/GF_SIN_BETA/coop_sqrt2       
       q = sqrt(3.d0-2.d0*GF_COS_BETA**2)
       g =  - sqrt(coop_pio2)/q**3*u*GF_COS_BETA*(1.d0+erf(t/q))*exp(-(u*GF_COS_BETA/q)**2) &
            - sqrt(coop_pio2)*u*GF_COS_BETA*(3.d0*GF_SIN_BETA**2+(u*GF_COS_BETA)**2-1.d0)*(1.d0+erf(t)) &
            - GF_SIN_BETA*exp(-t**2)*(1.d0/q**2 - 1 + 2.d0*GF_SIN_BETA**2 + (u*GF_COS_BETA)**2)
    case(3)
       t = -u*GF_COS_BETA
       q = 5.d0*GF_SIN_BETA
       if(t + q .lt. 0.d0 .and. t .lt. -3.d0)then
          g = -sqrt(coop_2pi)*t**2*(t**2+(6.d0*GF_SIN_BETA**2-3.d0))
       else
          g = coop_integrate(dgdv, min(t, 0.d0)-q, min(0.d0, t+q), 1.d-7) / GF_SIN_BETA
       endif
    case default
       stop "gaussian_peak_G only support dim = 1, 2, 3"
    end select
    g= g*exp(-u**2/2.d0)
  contains
    function dGdv(v)
      COOP_REAL dGdv, v
      dGdv =  coop_gaussian_peak_f(v, args) * v * exp(- (((v + u*GF_COS_BETA)/GF_SIN_BETA)**2)/2.d0)
    end function dGdv
  end function coop_gaussian_peak_Gv


  function coop_gaussian_peak_Gv_slow(u, args) result(g)
    COOP_REAL t, q, g, u
    type(coop_arguments)::args    
    t = -u*GF_COS_BETA
    q = 6.d0*GF_SIN_BETA
    g = coop_integrate(dgdv, min(t, 0.d0)-q, min(0.d0, t+q), 1.d-8) / GF_SIN_BETA*exp(-u**2/2.d0)
  contains
    function dGdv(v)
      COOP_REAL dGdv, v
      dGdv =  coop_gaussian_peak_f(v, args) * v * exp(- (((v + u*GF_COS_BETA)/GF_SIN_BETA)**2)/2.d0)
    end function dGdv
  end function coop_gaussian_peak_Gv_slow  
  

  function coop_gaussian_peak_Gu(u, args) result(g)
    COOP_REAL u, g
    type(coop_arguments)::args
    g= u * coop_gaussian_peak_G(u, args)
  end function coop_gaussian_peak_Gu  
  

  !!before marginalization of e
  function coop_gaussian_npeak_differential_2D_full(u, v, e, args) result(dnpk)
    COOP_REAL u, v, e, dnpk
    type(coop_arguments)::args
    dnpk = (GF_CONST*2.d0) * v**4*e*(1.d0-e**2)*exp(-(v*e)**2 - (u**2 + ((v + u*GF_COS_BETA)/GF_SIN_BETA)**2)/2.d0) 
  end function coop_gaussian_npeak_differential_2D_full

  function coop_gaussian_nmax(nu, args) result(nmax)
    COOP_REAL nmax, nu
    type(coop_arguments)::args
    if(nu .gt. 3.d0)then
       nmax = GF_CONST*coop_integrate(coop_gaussian_peak_G, nu, nu + 3.d0, args, 1.d-7)       
    else
       nmax = GF_CONST*coop_integrate(coop_gaussian_peak_G, max(nu, -5.d0), 6.d0, args, 1.d-7)
    endif
  end function coop_gaussian_nmax

  function coop_gaussian_peak_intu(nu, args) result(intu)
    COOP_REAL intu, nu
    type(coop_arguments)::args
    if(nu .gt. 3.d0)then
       intu = GF_CONST*coop_integrate(coop_gaussian_peak_Gu, nu, nu + 3.d0, args, 1.d-7)       
    else
       intu = GF_CONST*coop_integrate(coop_gaussian_peak_Gu, max(nu, -5.d0), 6.d0, args, 1.d-7)
    endif

  end function coop_gaussian_peak_intu

  function coop_gaussian_peak_intv(nu, args) result(intv)
    COOP_REAL intv, nu
    type(coop_arguments)::args
    if(nu .gt. 3.d0)then
       intv = GF_CONST*coop_integrate(coop_gaussian_peak_Gv, nu, nu + 3.d0, args, 1.d-7)       
    else
       intv = GF_CONST*coop_integrate(coop_gaussian_peak_Gv, max(nu, -5.d0), 6.d0, args, 1.d-7)
    endif
  end function coop_gaussian_peak_intv
  
  subroutine coop_gaussian_get_nonoriented_stacking_weights(nu, args, w0, w2)
    type(coop_arguments)::args
    COOP_REAL nu
    COOP_REAL:: A(2,2), nmax, vec(2), w0, w2
    A(1,1) = 1.d0
    A(2,2) = 1.d0
    A(1,2) = GF_COS_BETA
    A(2,1) = GF_COS_BETA
    nmax = coop_gaussian_nmax(nu, args)
    vec(1) = coop_gaussian_peak_intu(nu, args)/nmax
    vec(2) = coop_gaussian_peak_intv(nu, args)/nmax
    vec = matmul(A, vec)/GF_SIN_BETA**2
    w0 = vec(1)
    w2 = vec(2)
  end subroutine coop_gaussian_get_nonoriented_stacking_weights
  
end module coop_gaussian_peak_stat_mod
