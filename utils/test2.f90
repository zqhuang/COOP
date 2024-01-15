module gutils
  use coop_wrapper_utils
  use coop_expint_mod  
  implicit none
#include "constants.h"

  COOP_REAL:: r
  COOP_REAL,parameter::C0 = 0.6303307007539d0
  COOP_INT,parameter::tab_n = 30000
  COOP_REAL::tab_lnx(tab_n), tab_lnbeta(tab_n), tab_lnbeta2(tab_n)
  COOP_REAL,  dimension(13), parameter:: coefs = (/ &
       -0.0005000886722098459d0, & 
       0.0035749585243020276d0, & 
       -0.011550555275443703d0, & 
       0.022388497813689085d0, & 
       -0.029222724394900604d0, & 
       0.027490121791118233d0, & 
       -0.0197349256846174d0, & 
       0.011640637732273086d0, & 
       -0.0064023539311030475d0, & 
       0.004075030604809154d0, & 
       -0.004030613541561846d0, & 
       0.009019150446682427d0, & 
       -0.13707783616386848d0 /)
  
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
    if(lnx .ge. tab_lnx(tab_n))then
       beta = exp(-x)*(x*(x+2.d0)-1.d0)/x**2/(x+3.d0)
    elseif(lnx .le. tab_lnx(1))then
       beta = 1.d0/x + 0.5d0*log(x) - 0.630330697184d0 - x/12.d0
    else
       call coop_splint(tab_n, tab_lnx, tab_lnbeta, tab_lnbeta2, lnx, beta)
       beta = exp(beta)
    endif
  end function beta_of_x

  

  function C_of_a(a)
    COOP_REAL::a, C_of_a, x
    C_of_a = C0 + diffC(a)
  end function C_of_a

  function diffC(a)  !!C(a) - C(0)
    COOP_REAL::a, diffC, a2
    COOP_INT::i
    a2 = a*a
    diffC =   a2*coefs(1)
    do i=2, size(coefs)
       diffC = a2*(diffC + coefs(i))
    enddo
  end function diffC
  
  function h1_of_t(t) result(h1)
    COOP_REAL::x, h1, t
    x = exp(t)*sqrt(r)
    if(x .gt. 1.d-6)then
       h1 = x*beta_of_x(x**2)-1.d0/x
    else  !!use more accurate asympototic
       h1 = x*(log(x) - C0 - x**2/12.d0)
    endif
  end function h1_of_t

  recursive function h2_of_t(t) result(h2)
    COOP_REAL::h2, t, tbyst
    if(t .lt. 0.d0) then
       h2 = h2_of_t(-t)
       return
    endif
    if(t .gt. 1.d-5)then
       tbyst = t/sinh(t)
    else
       tbyst = 1.d0/(1.d0+t**2/6.d0)
    endif
    h2 = (t+C0)*exp(-t) - diffC(exp(-2*t))*exp(t) + 1.5d0*tbyst !(2*t+ C0)*exp(-t) - diffC(exp(-2*t)) *exp(t) + tbyst*(1.d0+exp(-2*t)/2.d0)
  end function h2_of_t


  function h_of_t(t) result(h)
    COOP_REAL::t, h
    h = (h1_of_t(t) + h1_of_t(-t))/sqrt(r) + h2_of_t(t) 
  end function h_of_t


  function lnh_of_t(t) result(lnh)
    COOP_REAL::t, lnh
    lnh = log(-h_of_t(t))
  end function lnh_of_t
  

  subroutine dump_h(r, tmax, nt, filename)
    COOP_UNKNOWN_STRING::filename
    COOP_INT::nt, i
    COOP_REAL::tmax, t, lnr, r, tbyst, hmain, h
    lnr = log(r)
    open(31, file=filename)
    do i=1, nt
       t = tmax*(i-1)/(nt-1.d0)
       h = h_of_t(t)
       write(31, "(4E18.9)") t,  h
    enddo
    close(31)
  end subroutine dump_h


  function inth2_analytic(z) result(f)
    COOP_COMPLEX::mu, z, f
    COOP_INT::i
    mu = z -0.5d0
    f = 3.d0*coop_pi**2/4.d0/cos(coop_pi*mu)**2 + (2.d0+8.d0*mu**2)/(1.d0-4.d0*mu**2)**2 - 0.5*C0/(mu**2-0.25d0)
    do i=1, size(coefs)
       f = f + coefs(size(coefs)-i+1)*(8*i-2)/(4*mu**2 - (4*i-1.d0)**2)
    enddo
  end function inth2_analytic

end module gutils

program test
  use gutils
  use coop_wrapper_utils
  use coop_expint_mod
  implicit none
#include "constants.h"
  COOP_COMPLEX:: z, ind, inth, theory
  COOP_REAL,dimension(10),parameter::zeta_zeros = (/ 14.134725142d0, 21.022039639d0, 25.010857580d0, 30.424876126d0, 32.935061588d0,     37.586178159d0,     40.918719012d0,      43.327073281d0,      48.005150881d0,    49.773832478d0 /)
  COOP_INT::i, n

  do i=0, 5000
     print*, i*0.01d0, h2_of_t(i*0.01d0)*sinh(i*0.01d0*0.5d0)
  enddo
  stop

  call load_beta()
  r = 1.d-20  
  call dump_h(r, 18.d0-log(r)/2.d0, 10000, "h_20.txt")
  stop
  open(11, file="log1.txt")
  do i= 200, 4600
     z = cmplx(0.75d0, i*0.005d0 )
     ind = 2.d0*z - 1.d0
     inth = coop_periodsym_int(f = h2_of_t, a = 30.d0, z= ind, accuracy = 1.d-8)
     theory = inth2_analytic(z)
   !  theory =  ((coop_Riemannzeta(z)*exp(coop_complex_log_gamma(z) - z*log(r))/z + coop_RiemannZeta(1.d0-z)*exp(coop_complex_log_gamma(1.d0-z)+(z-1.d0)*log(r))/(1.d0-z)) + coop_RiemannZeta(z)*coop_RiemannZeta(1.d0-z)* coop_pi/ sin(coop_pi*z))/2.d0
     write(11, "(5E16.7)") aimag(z), inth , theory
     write(*, "(5E14.5)") aimag(z), inth , theory
  enddo
  close(11)
  
end program test
