program test
  use coop_wrapper_utils
  use coop_expint_mod
  implicit none
#include "constants.h"
  character(LEN=*),parameter::postfix = "_h20_s2.txt"
  COOP_REAL::r
  COOP_INT,parameter::ntab = 10000
  COOP_COMPLEX:: z, ind, inth
  COOP_REAL,dimension(10),parameter::zeta_zeros = (/ 14.134725142d0, 21.022039639d0, 25.010857580d0, 30.424876126d0, 32.935061588d0,     37.586178159d0,     40.918719012d0,      43.327073281d0,      48.005150881d0,    49.773832478d0 /)
  COOP_REAL::tab_lnl(ntab), tab_lnh(ntab), tab_lnh2(ntab), tmp1, tmp2, dlnl
  COOP_REAL::left_slope, right_slope
  COOP_INT::i, n
  open(11, file="tab"//postfix)
  read(11, *) r
  do i = 1, ntab
     read(11, *) tmp1, tmp2, tab_lnl(i), tab_lnh(i)
  enddo
  close(11)

  
  
  call coop_spline(ntab, tab_lnl, tab_lnh, tab_lnh2, -1.5d0, 0.d0)
  open(11, file="log"//postfix)
  do i= 2200, 3400
     z = cmplx(0.5d0, i*0.005d0 )
     ind = 2.d0*z - 1.d0
     inth = coop_period_int(f = h_of_t, a = -60.d0, b = 60.d0, z= ind, accuracy = 1.d-9)
     write(11, "(5E16.7)") aimag(z), inth, ((coop_Riemannzeta(z)*exp(coop_complex_log_gamma(z) - z*log(r))/z + coop_RiemannZeta(1.d0-z)*exp(coop_complex_log_gamma(1.d0-z)+(z-1.d0)*log(r))/(1.d0-z)) + coop_RiemannZeta(z)*coop_RiemannZeta(1.d0-z)* coop_pi/ sin(coop_pi*z))/2.d0
     if(mod(i, 100).eq.0)write(*, "(5E16.7)") aimag(z), inth, ((coop_Riemannzeta(z)*exp(coop_complex_log_gamma(z) - z*log(r))/z + coop_RiemannZeta(1.d0-z)*exp(coop_complex_log_gamma(1.d0-z)+(z-1.d0)*log(r))/(1.d0-z)) + coop_RiemannZeta(z)*coop_RiemannZeta(1.d0-z)* coop_pi/ sin(coop_pi*z))/2.d0     
  enddo
  close(11)
contains

  function h_of_t(t) result(h)
    COOP_REAL::lnl, rl, h, t
    COOP_INT::il
    lnl = -2.d0*abs(t) - log(r)
    if(lnl .lt. tab_lnl(1))then
       h = tab_lnh(1) + 1.5*(lnl - tab_lnl(1))
    else
       call coop_splint(ntab, tab_lnl, tab_lnh, tab_lnh2, lnl, h)
    endif
    h = - exp(h) /sqrt(r)
  end function h_of_t

  function gauss(t)
    COOP_REAL::t, gauss
    gauss = exp(-t**2/2.d0)
  end function gauss
  
end program test
