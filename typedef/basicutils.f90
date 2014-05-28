module coop_basicutils
  use coop_constants
  use coop_svd
  implicit none
#include "constants.h"

contains
  
  subroutine coop_set_uniform(n, x, lower, upper, logscale)
    COOP_INT n, i
    COOP_REAL x(n), lower, upper, rlow, dx
    logical,optional::logscale
    x(1) = lower
    x(n) = upper
    if(present(logscale))then
       if(logscale)then
          dx = (log(upper) - log(lower))/(n-1)
          rlow = log(lower) - dx
          !$omp parallel do
          do i = 2, n-1
             x(i) = exp(rlow + dx*i)
          enddo
          !$omp end parallel do
          return
       endif
    endif
    dx = (upper-lower)/(n-1)
    rlow = lower-dx
    !$omp parallel do
    do i = 2, n-1
       x(i) =rlow + dx*i
    enddo
    !$omp end parallel do

  end subroutine coop_set_uniform

  subroutine coop_locate(n, x, needle, loc, res)
    COOP_INT n, loc, imin, imax
    COOP_REAL x(n), needle, res
    if(x(1).le. x(n))then
       if(needle .lt. x(1))then
          loc = 0
          return
       endif
       if(needle .gt. x(n))then
          loc = n
          return
       endif
       imin = 1
       imax = n
       do while(imax - imin .gt. 1)
          loc = (imax + imin)/2
          if(x(loc).le. needle)then
             imin = loc
          else
             imax = loc
          endif
       enddo
       loc = imin
       if(needle .gt. x(imin))then
          res = (needle - x(imin))/(x(imax)-x(imin))
       else
          res = 0
       endif
       return
    else
       if(needle .gt. x(1))then
          loc = 0
          return
       endif
       if(needle .lt. x(n))then
          loc = n
          return
       endif
       imin = 1
       imax = n
       do while(imax - imin .gt. 1)
          loc = (imax + imin)/2
          if(x(loc).ge. needle)then
             imin = loc
          else
             imax = loc
          endif
       enddo
       loc = imin
       if(needle .lt. x(imin))then
          res = (needle - x(imin))/(x(imax)-x(imin))
       else
          res = 0
       endif
       return       
    endif
  end subroutine coop_locate


  subroutine coop_spline(n, x, y, y2, ypl, ypr)
    COOP_INT n, i
    COOP_REAL x(n), y(n), y2(n)
    COOP_REAL, optional::ypl,ypr
    COOP_REAL yil, yir, bet, dxr, dxl
    COOP_REAL gam(n-1)
    if(n.le.1)then
       y2 = 0.
       return
    endif
    dxr = x(2) - x(1)
    yir=(y(2)-y(1))/dxr
    if(present(ypl))then
       y2(1)=(yir-ypl)/dxr*3.
       gam(1)= 0.5
    else
       y2(1)=0.
       gam(1)=0.
    endif
    dxr = dxr/6.
    do i=2, n-1
       dxl = dxr
       dxr=x(i+1)-x(i)
       bet=(x(i+1)-x(i-1))/3.-dxl*gam(i-1)
       if(abs(bet) .lt. 1.d-30) stop 'Error in SPLinE.'
       yil=yir
       yir=(y(i+1)-y(i))/dxr
       y2(i)=(yir-yil-dxl*y2(i-1))/bet
       dxr=dxr/6.
       gam(i)=dxr/bet
    enddo
    if(present(ypr))then
       bet=(x(n)-x(n-1))/3.-dxr*gam(n-1)
       if(abs(bet) .lt. 1.d-30) stop 'Error in SPLinE.'
       y2(n)=(ypr-yir-dxr*y2(n-1))/bet
    else
       y2(n)=0.
    endif
    do i=n-1, 1 , -1
       y2(i)=y2(i)-gam(i)*y2(i+1)
    enddo
  end subroutine coop_spline


  subroutine coop_splint(n, x, y, y2, xs, ys)
    COOP_INT n, l, r
    COOP_REAL x(n), y(n), y2(n)
    COOP_REAL xs, ys, a, b
    call coop_locate(n, x, xs, l, b)
    if(l .lt. 1)then
       ys=y(1)
       return
    endif
    if( l .ge. n)then
       ys=y(n)
       return
    endif
    r = l + 1
    a =  1. - b
    ys=y(l)*a+y(r)*b+  &
         (y2(l)*(a*a-1.)*a+y2(r)*(b*b-1.)*b)/6.*(x(r)-x(l))**2
  end subroutine coop_splint


  subroutine coop_spline_uniform(n, y, y2, ypl, ypr)
    COOP_INT n, i
    COOP_REAL  y(n), y2(n)
    COOP_REAL, optional::ypl,ypr
    COOP_REAL yil, yir, bet, dxl, dxr
    COOP_REAL gam(n-1)
    if(n.le.1)then
       y2 = 0.
       return
    endif
    dxr = 1.
    yir=y(2)-y(1)
    if(present(ypl))then
       y2(1)=(yir-ypl)*3.
       gam(1)= 0.5
    else
       y2(1)=0.
       gam(1)=0.
    endif
    dxr = dxr/6.
    do i=2, n-1
       dxl = dxr
       dxr=1.
       bet=1.d0/3.-dxl*gam(i-1)
       if(abs(bet) .lt. 1.d-30) stop 'Error in SPLinE.'
       yil=yir
       yir=(y(i+1)-y(i))/dxr
       y2(i)=(yir-yil-dxl*y2(i-1))/bet
       dxr=dxr/6.
       gam(i)=dxr/bet
    enddo
    if(present(ypr))then
       bet=1.d0/3.-dxr*gam(n-1)
       if(abs(bet) .lt. 1.d-30) stop 'Error in SPLinE.'
       y2(n)=(ypr-yir-dxr*y2(n-1))/bet
    else
       y2(n)=0.
    endif
    do i=n-1, 1 , -1
       y2(i)=y2(i)-gam(i)*y2(i+1)
    enddo
    y2 =  y2/6.
  end subroutine coop_spline_uniform

  subroutine coop_splint_uniform(n, xmin, dx, y, y2, xs, ys)
    COOP_INT n, l, r
    COOP_REAL xmin, dx, y(n), y2(n), xs, ys, a, b
    b = (xs - xmin)/dx + 1.
    l = floor(b)
    if(l.lt.1)then
       ys = y(1)
       return
    endif
    if(l.ge.n)then
       ys = y(n)
       return
    endif
    b=b-l
    r = l+1
    a=1.-b
    ys=y(l)*a+y(r)*b+  &
         (y2(l)*(a*a-1.)*a+y2(r)*(b*b-1.)*b)
  end subroutine coop_splint_uniform

  subroutine coop_cheb_eval_all(n, x, y)
    COOP_INT,intent(IN):: n
    COOP_REAL,intent(IN):: x
    COOP_REAL,intent(OUT)::y(0:n)
    COOP_REAL twox
    COOP_INT i
    y(0) = 1.d0
    if(n.eq.0) return
    y(1) = x
    twox = 2.*x
    do i = 2, n
       y(i) = twox * y(i-1) - y(i-2)
    enddo
  end subroutine coop_cheb_eval_all


  subroutine coop_get_cheb_value(n, c, x, y)
    COOP_INT n
    COOP_REAL c(n), x, y
    COOP_REAL twox, y1, y2, y3
    COOP_INT i
    twox = x+x
    y3 = 0.d0
    y2 = 0.d0
    do i = n, 2, -1
       y1 = c(i) + twox * y2 - y3
       y3 = y2
       y2 = y1
    enddo
    y = c(1) + x* y2 - y3    
  end subroutine coop_get_cheb_value


  subroutine coop_fit_template(n, m, y, tpls, c)
    COOP_INT,intent(IN):: n, m
    COOP_REAL,intent(IN)::y(n)
    COOP_REAL,intent(INOUT)::tpls(n, m)
    COOP_REAL,intent(OUT):: c(m)
    call coop_svd_least_square_one(n, m, tpls, y, c)
  end subroutine coop_fit_template


  subroutine coop_chebfit(n, x,y, m, a,b,c)
    COOP_INT,intent(in):: n, m
    COOP_REAL,intent(in)::a, b
    COOP_REAL,intent(in)::x(n), y(n)
    COOP_REAL,intent(out)::c(m)
    COOP_REAL fx(n, m ),t(n)
    COOP_INT i
    if(m.gt.n)then
       write(*,*) "coop_chebfit: Not enough data for chebyshev polynomial fit."
       stop
    endif
    t=2.d0*(x-a)/(b-a)-1.d0
    !$omp parallel do
    do i=1,n
       call coop_cheb_eval_all(m-1,t(i), Fx(i,1:m))
    enddo
    !$omp end parallel do
    call coop_fit_template(n, m, y, fx, c)
  end subroutine coop_chebfit

  subroutine coop_chebfit_uniform(n, y, m, c)
    COOP_INT,intent(in):: n, m
    COOP_REAL,intent(in)::y(n)
    COOP_REAL,intent(out)::c(m)
    COOP_REAL fx(n, m ),dt
    COOP_INT i
    if(m.gt.n)then
       write(*,*) "coop_chebfit_uniform: Not enough data for chebyshev polynomial fit."
       stop
    endif
    dt = 2.d0/(n-1)
    !$omp parallel do
    do i=1,n
       call coop_cheb_eval_all(m-1, -1.+ (i-1) * dt, Fx(i,1:m))
    enddo
    !$omp end parallel do
    call coop_fit_template(n, m, y, fx, c)
  end subroutine coop_chebfit_uniform

  subroutine coop_chebeval(n, a, b, c, x, y)
    COOP_INT n
    COOP_REAL a, b, c(n), x, y, t
    t = 2.*(x-a)/(b-a)-1.
    call  coop_get_cheb_value(n, c, t, y)
  end subroutine coop_chebeval

end module coop_basicutils
