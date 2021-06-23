module tmp
#include "constants.h"      
  use coop_wrapper_utils

contains

  subroutine fcn(n, t, y, yp)
    COOP_INT::n
    COOP_REAL::t, y(n), yp(n)
    COOP_REAL,parameter::k = 1.d2
    COOP_REAL,parameter::m = 0.27
    COOP_REAL,parameter::g = 10.d0
    COOP_REAL,parameter::y0 = (m*g/k)**(1.d0/3.d0)
    yp(1) = y(2)
    yp(2) = -k/m*(y(1)+y0)**3 + g
  end subroutine fcn

  function kessence_F(x) result(F)
    real*8::x, F, y
    y = x**3
    if(x > 0.03d0)then
       F = sqrt(1.d0+1.d0/y) - log(sqrt(y) + sqrt(1.d0+y))/y
    else
       F = ((2.d0/3.d0) + y*(- 0.2d0+y*(3.d0/28.d0)))*sqrt(y)
    endif
  end function kessence_F

  function kessence_dFdx(x) result(dFdx)
    real*8::x, dFdx, y,sqrty, sqrtz
    y = x**3
    if(x > 0.03d0)then
       sqrty = sqrt(y)
       sqrtz = sqrt(1.d0+y)
       dFdx = 3.d0/x/sqrty*(log(sqrty+sqrtz)/sqrty-1.d0/sqrtz)
    else
       dFdx = (1.d0 + y*( - 0.9d0 + y* (45.d0/56.d0)))*sqrt(x)
    endif
  end function kessence_dFdx


  function wrat(x)
    real*8 x, wrat, y
    if(x > 0.3d0)then
       wrat = -2*x*kessence_dFdx(x)/kessence_F(x)
    else
       y = x**3
       wrat = -3.d0+y*(1.8d0+y*(-243.d0/175.d0 + y * (1023.d0/875.d0- y*(346431.d0/336875.d0))))
    end if
  end function wrat


  function inverse_wrat(wr) result(x)
    real*8 wr, x, c
    real*8,parameter::a = 243.d0/175.d0, b=1.8d0
    real*8,parameter::xpiv=1.4d0, wrpiv = -1.3155941865393341, slop_piv = 1.34376
    c = wr + 3.d0
    if(c < 1.d-5)then
       x = (c/b)**(1.d0/3.d0)
    else if(c < 0.5d0)then
       x = ((b - sqrt(b**2-4*a*c))/(2*a))**(1.d0/3.d0)
    else
       x = (wr-wrpiv)/slop_piv + xpiv
    endif
  end function inverse_wrat
    

end module tmp



program Test
#include "constants.h"    
  use coop_wrapper_utils
  use tmp
  implicit none
  integer,parameter::n=201
  integer i, j
  real*8 x(n), wr(n), xinv(n)
  type(coop_asy)::fig
  call fig%open("wrat.txt")
  call fig%init(xlabel="$x$", ylabel="$w$")
  call coop_set_uniform(n, x, 0.d0, 2.5d0)
  do i=1, n
     wr(i) = wrat(x(i))
     xinv(i) = inverse_wrat(wr(i))
  enddo
  call fig%plot(x, xinv-x)
  call fig%close()
end program Test



