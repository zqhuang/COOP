module gutils
  use coop_wrapper_utils
  implicit none
#include "constants.h"

  COOP_INT,parameter::ntab = 103308
  COOP_REAL::tab_lnx_start, tab_dlnx, tab_lnx_end
  COOP_REAL,dimension(ntab)::tab_lng0, tab_lng1
  COOP_REAL,parameter:: r = 1.d-20
  COOP_REAL,parameter:: sqrtr = sqrt(r)
  
contains

  function g0_small(lnx) result(g0)
    COOP_REAL::lnx, x, g0
    x = dexp(lnx)
    g0 = 1.d0 + x*(lnx/2.d0 - x/12.d0 - 0.630330697184d0)
  end function g0_small

  function g0_large(lnx) result(g0)
    COOP_REAL::x, g0, lnx
    x = dexp(lnx)
    g0 = (1.d0-(1.d0-(2.d0-6.d0/x)/x)/x)*dexp(-x) 
  end function g0_large

  function g1_small(lnx) result(g1)
    COOP_REAL::lnx, x, g1
    x = dexp(lnx)
    g1 = -lnx + x*(0.5d0 - x/24.d0)
  end function g1_small

  function g1_large(lnx) result(g1)
    COOP_REAL::ex, g1, lnx
    ex = dexp(-dexp(lnx))
    g1 = ex*(1.d0+ex*(0.5d0 + ex*(1.d0/3.d0 + ex/4.d0)))
  end function g1_large

  function g0_of_lnx(lnx) result(g0)
    COOP_REAL::lnx, g0, ir
    COOP_INT::i
    if(lnx .le. tab_lnx_start)then
       g0 = g0_small(lnx)
    elseif( lnx .ge. tab_lnx_end)then
       g0 = g0_large(lnx)
    else
       ir = (lnx - tab_lnx_start)/tab_dlnx+1.d0
       i = floor(ir)
       ir = ir - i
       g0 = dexp(tab_lng0(i)*(1.d0-ir) + tab_lng0(i+1)*ir)
    endif
  end function g0_of_lnx

  function g1_of_lnx(lnx) result(g1)
    COOP_REAL::lnx, g1, ir
    COOP_INT::i
    if(lnx .le. tab_lnx_start)then
       g1 = g1_small(lnx)
    elseif( lnx .ge. tab_lnx_end)then
       g1 = g1_large(lnx)
    else
       ir = (lnx - tab_lnx_start)/tab_dlnx+1.d0
       i = floor(ir)
       ir = ir - i
       g1 = dexp(tab_lng1(i)*(1.d0-ir) + tab_lng1(i+1)*ir)
    endif
  end function g1_of_lnx

  subroutine load_gs()
    COOP_INT::i
    COOP_REAL::x
    open(11, file="g0.txt")
    do i = ntab, 1, -1
       read(11, *) x, tab_lng0(i)
    enddo
    close(11)
    open(11, file="g1.txt")
    do i = ntab, 1, -1
       read(11, *) x, tab_lng1(i)
    enddo
    close(11)
    tab_dlnx = -dlog(1.d0-1.d0/8192.d0)
    tab_lnx_end = dlog(30.d0)
    tab_lnx_start =  tab_lnx_end - tab_dlnx *(ntab-1)
    tab_lnx_end = tab_lnx_end - 1.d-8
  end subroutine load_gs

  function zeta_int(t) result(f)
    COOP_REAL::t, x, f, tt, tt2, lnx
    tt = dexp(t)
    tt2 = tt**2
    x = r/tt2
    lnx = dlog(x)
    f = (g0_of_lnx(lnx)-1.d0)*tt/sqrtr
  end function zeta_int

  recursive function zeta_int2(t) result(f) !!tan theta  = exp(t)
    COOP_REAL::t, x, f, tt, tt2, lnx
    if(t > 0.d0)then
       f = zeta_int2(-t)
       return
    endif
    tt = dexp(t)
    tt2 = tt**2
    x = r/tt2
    lnx = dlog(x)
    if(tt2 < 0.9999d0)then
       f = (g0_of_lnx(lnx)-1.d0 + r/2.d0 * (-g1_of_lnx(lnx)-dlog(tt2)*(2.d0+tt2)/(1.d0-tt2)))*tt/sqrtr
    else
       f = (g0_of_lnx(lnx)-1.d0 + r/2.d0 * (-g1_of_lnx(lnx)+(2.d0+tt2)*(1.d0+(1.d0-tt2)*(0.5d0+(1-tt2)/3.d0))))*tt/sqrtr
    endif
  end function zeta_int2
  
end module gutils

program test
  use gutils
  use coop_expint_mod
  use coop_wrapper_utils
  implicit none
  COOP_REAL,dimension(10),parameter::zeta_zeros = (/ 14.134725142d0, 21.022039639d0, 25.010857580d0, 30.424876126d0, 32.935061588d0,     37.586178159d0,     40.918719012d0,      43.327073281d0,      48.005150881d0,    49.773832478d0 /)
  COOP_COMPLEX::z, zint
  COOP_INT::i
  COOP_REAL::t
  call load_gs()
  do i=1000, 0, -1
     t = i*0.03d0
     print*, -2*t-log(r),log(- zeta_int(-t)*sqrtr)-t
  enddo
end program test
