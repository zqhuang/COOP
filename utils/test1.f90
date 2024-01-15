module gutils
  use coop_wrapper_utils
  use coop_expint_mod  
  implicit none
#include "constants.h"

  COOP_REAL:: param_lam !!lam = e^{-2t}
  COOP_INT,parameter::tab_n = 30000
  COOP_REAL::tab_lnx(tab_n), tab_lnbeta(tab_n), tab_lnbeta2(tab_n)
  
contains

  subroutine load_beta()
    COOP_INT::i
    open(21, file="logbeta_logx.txt")
    do i=1, tab_n
       read(21, *) tab_lnx(i), tab_lnbeta(i)
    enddo
    close(21)
    call coop_spline(tab_n, tab_lnx, tab_lnbeta, tab_lnbeta2)
  end subroutine load_beta

  function beta_of_x(x) result(beta)
    COOP_REAL::x, beta, lnx
    lnx = log(x)
    if(x > 50.d0)then
       beta = exp(-x)*(x*(x+2.d0)-1.d0)/x**2/(x+3.d0)
    elseif(lnx .le. tab_lnx(1))then
       beta = 1.d0/x + 0.5d0*log(x) - 0.630330697184d0 - x/12.d0
    else
       call coop_splint(tab_n, tab_lnx, tab_lnbeta, tab_lnbeta2, lnx, beta)
       beta = exp(beta)
    endif
  end function beta_of_x
  
  function h1_integrand(s) result(h1)
    COOP_REAL::h1, s, st, s_err, st_err
    st = s*param_lam
    if(s .lt. 5.d-3)then
       s_err = 0.5d0 + s*(1.d0/6.d0 + s * (1.d0/24.d0 + s * (1.d0/120.d0 + s / 720.d0)))
       st_err = 0.5d0 + st*(1.d0/6.d0 + st * (1.d0/24.d0 + st * (1.d0/120.d0 + st / 720.d0)))
       h1 = (s_err + (param_lam + st*s_err)*st_err)/(1.d0 + s*s_err)/(1.d0+st*st_err)/st
    else if(st .lt. 5.d-3)then
       st_err = 0.5d0 + st*(1.d0/6.d0 + st * (1.d0/24.d0 + st * (1.d0/120.d0 + st / 720.d0)))
       h1 = (1.d0/s - 1.d0/(exp(s)-1.d0)/(1.d0+st*st_err))/st
    else
       h1 = 1.d0/s/st - 1.d0/(exp(s)-1.d0)/(exp(st)-1.d0)
    endif
  end function h1_integrand


  function h2_integrand(s) result(h2)
    COOP_REAL::h2, s, st, s_err, st_err, b
    st = s*(1.d0-param_lam)
    if(s .lt. 5.d-3)then
       s_err = 0.5d0 + s*(1.d0/6.d0 + s * (1.d0/24.d0 + s * (1.d0/120.d0 + s / 720.d0)))
       st_err = 0.5d0 + st*(1.d0/6.d0 + st * (1.d0/24.d0 + st * (1.d0/120.d0 + st / 720.d0)))
       b = s*param_lam
       b = param_lam*(1.d0-b*(0.5d0 - b*(1.d0/6.d0 - b*(1.d0/24.d0 - b/120.d0))))
       h2 =(s_err + (1.d0-param_lam + st*s_err)*st_err + b)/(1.d0 + s*s_err)/(1.d0+st*st_err)/st
    else if(st .lt. 5.d-3)then
       st_err = 0.5d0 + st*(1.d0/6.d0 + st * (1.d0/24.d0 + st * (1.d0/120.d0 + st / 720.d0)))
       h2 = (1.d0/s - exp(-s*param_lam)/(exp(s)-1.d0)/(1.d0+st*st_err))/st
    else
       h2 = 1.d0/s/st - exp(-s)/(exp(s)-1.d0)/(1.d0-exp(-st))
    endif
  end function h2_integrand


  function h3_integrand(s) result(h3)
    COOP_REAL::h3, s, st, s_err, st_err, b, ss_err
    st = s*(1.d0-param_lam)
    ss_err = s*(1.d0/6.d0 + s * (1.d0/24.d0 + s * (1.d0/120.d0 + s / 720.d0)))
    s_err = 0.5d0 + ss_err
    st_err = 0.5d0 + st*(1.d0/6.d0 + st * (1.d0/24.d0 + st * (1.d0/120.d0 + st / 720.d0)))
    b = s*param_lam
    b = param_lam*(1.d0-b*(0.5d0 - b*(1.d0/6.d0 - b*(1.d0/24.d0 - b/120.d0))))
    h3 =(ss_err - 0.5d0*(s*s_err + st*st_err*(1.+s*s_err)) + (1.d0-param_lam + st*s_err)*st_err + b)/(1.d0 + s*s_err)/(1.d0+st*st_err)/s
  end function h3_integrand
  

  function h_func(r) result(h)
    COOP_REAL::r, h
    COOP_REAL,parameter::accuracy= 1.d-6
    COOP_REAL,parameter::upper = 20.d0
    COOP_REAL::corr
    corr = 1.d0/(r + upper*param_lam)
    if(param_lam .lt. 0.2d0)then
       h = sqrt(param_lam)*(coop_expint(h2_integrand, r, r/param_lam, accuracy) - coop_expint(h1_integrand, r/param_lam, r/param_lam + upper, accuracy)-corr)
    else if(param_lam .lt. 0.999d0)then
       h = sqrt(param_lam)*(coop_integrate(h2_integrand, r, r/param_lam, accuracy) - coop_expint(h1_integrand, r/param_lam, r/param_lam + upper, accuracy)-corr)       
    else
       h = sqrt(param_lam)*(coop_integrate(h3_integrand, r, r/param_lam, accuracy)/(1.d0-param_lam)+0.5d0*(1.d0+(1.d0-param_lam)/2.d0) - coop_expint(h1_integrand, r/param_lam, r/param_lam + upper, accuracy)-corr)       
    endif
  end function h_func


  function c1int(s)
    COOP_REAL::s, c1int, err, tmp
    if(s .lt. 1.d-2)then
       tmp = 1.d0/6.d0+s*(1.d0/24.d0+s*(1.d0/120.d0+s*(1.d0/720.d0+s/5040.d0)))
       err = 0.5d0+tmp*s
       c1int = (2*err - tmp*2.d0 - err**2*(1.d0-s))/(1.d0+err*s)**2
    else
       c1int = -1.d0/s**2 + 1.d0/s + 1.d0/(exp(s)-1.d0)**2
    endif
  end function c1int

  function c0int(s)
    COOP_REAL::s, c0int, err, tmp
    if(s .lt. 1.d-2)then
       tmp = 1.d0/6.d0+s*(1.d0/24.d0+s*(1.d0/120.d0+s*(1.d0/720.d0+s/5040.d0)))
       err = 0.5d0+tmp*s
       c0int = (err/2.d0-tmp)/(1.d0+err*s)
    else
       c0int = (-1.d0/s + 0.5d0 + 1.d0/(exp(s)-1.d0))/s
    endif
  end function c0int


  function caint(s)
    COOP_REAL,parameter::a = 0.1d0
    COOP_REAL::s, caint, err1, tmp1, err2, tmp2, as
    if(s .lt. 1.d-2)then
       tmp1 = 1.d0/6.d0+s*(1.d0/24.d0+s*(1.d0/120.d0+s*(1.d0/720.d0+s/5040.d0)))
       err1 = 0.5d0+tmp1*s
       as = a*s
       tmp2 = (1.d0/6.d0+as*(1.d0/24.d0+as*(1.d0/120.d0+as*(1.d0/720.d0+as/5040.d0))))*a**2
       err2 = 0.5d0*a + tmp2*s
       caint = ((1.d0+a)/2.d0*(err1 + err2*(1.d0+err1*s))- tmp1 - tmp2 - err1*err2)/(1.d0+err1*s)/(1.d0+err2*s)
    else
       if(a .eq. 0.d0)then
          caint = (-1.d0/s + (1.d0+a)/2.d0)/s + 1.d0/(exp(s)-1.d0)/s
       else
          caint = (-1.d0/s + (1.d0+a)/2.d0)/s + a/(exp(s)-1.d0)/(exp(a*s)-1.d0)        
       endif
    endif
  end function caint

  function C_of_a(a)
    COOP_REAL::a, C_of_a
    C_of_a =0.630330697 + a* (-0.00000489 + a* (-0.13698186+ a*(-0.00070399+a*(0.011603725 + a*(-0.005127467 + 0.000883786*a)))))
  end function C_of_a

  function h_approx(r)
    COOP_REAL::t, h_approx, tbyst, r
    t = -log(param_lam)/2.d0
    if(r/param_lam .lt. 1.d-3)then
       if(t>1.d-5)then
          tbyst = t/sinh(t)
       else
          tbyst = 1.d0/(1.d0+t**2/6.d0)
       endif
       h_approx = (2.d0*t+log(r))*cosh(t) - C_of_a(param_lam)*exp(t) + tbyst*(1.d0+param_lam/2.d0) - (3.d0*exp(t)+5.d0*exp(-t)+exp(-3.d0*t))/6.d0*r
    else
       h_approx = exp(t)*beta_of_x(r/param_lam) + (0.5d0 - log(r) - 1.d0/r)*exp(-t)
    endif
  end function h_approx
  
  

end module gutils

program test
  use gutils
  use coop_wrapper_utils
  implicit none
  COOP_INT,parameter::n = 30
  COOP_REAL::lam(n), r, h, happ, t
  COOP_INT::i
  COOP_REAL,parameter::upper = 20.d0
  r = 1.d-20
  call load_beta()
!  do i=0, 10
!     a = 0.9 + 0.01d0 * i
!     write(*, "(2F18.11)") a,  -(coop_expint(caint, r, upper, 1.d-7) -1.d0/upper - (1.d0+a)/2.d0*log(upper) + (a**2+3*a+1)/12.d0*r) ! - (4.e-8+3.6e-8*a)
!  enddo
!  stop
  call coop_set_uniform(n, lam, r, r*1.d4, logscale = .true.)
!  open(11, file="tab_r"//COOP_STR_OF(log10(r))//".txt")
!  write(11, *) r
  do i=1, n
     param_lam = lam(i)
     t = -log(param_lam)/2.d0
     h = h_func(r)
     !     happ = (2.d0*t+log(r))*cosh(t) - (0.63-0.13*param_lam**2)*exp(t)+t/sinh(t)*(1.d0+param_lam/2.d0)
     happ = h_approx(r)
     print*, r/lam(i), h, happ, beta_of_x(r/param_lam), beta_of_x(r/param_lam)/sqrt(param_lam), (0.5d0 - log(r) - 1.d0/r)*sqrt(param_lam) !log(lam(i)/r), log(-h*sqrt(r))
     !write(11, *) lam(i), h, log(lam(i)/r), log(-h*sqrt(r))
  enddo
!  close(11)
end program test
