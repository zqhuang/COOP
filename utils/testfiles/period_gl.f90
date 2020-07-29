

module tmp
  use coop_wrapper_utils
#include "constants.h"
  real*8,parameter::theta0 = coop_pio2
contains


  function sinc(x)
    COOP_REAL::sinc, x, x2
    if(abs(x).gt. 1.d-3)then
       sinc = sin(x)/x
    else
       x2 = x*x
       sinc = 1.d0+x2*(-1.d0/6.d0+ x2*(1.d0/120.d0+x2*(-1.d0/5040.d0+x2*(1.d0/(72.d0*5040.d0)))))
    endif
  end function sinc
  
  function intdt(theta)
    real*8 theta, intdt
    intdt = 1.d0/sqrt(sin((theta0-theta)/2.d0) * sin((theta0+theta)/2.d0))
  end function intdt

  function intdx(x)
    COOP_REAL::x,intdx
    intdx = 1.d0/sqrt(sinc(x**2)*sin(theta0-x**2))
  end function intdx
  
  function simple_integrate(f, a, b) result(s)
    COOP_INT,parameter::n = 1000000
    COOP_REAL::f, a, b, s1,s2, dx, x, s
    external f
    COOP_INT::i
    s1 = 0.d0
    s2 = 0.d0
    dx = (b-a)/n
    do i=1, n-1
       s1 = s1 + f(a + dx*i)
       s2 = s2 + f(a + dx*(i-0.5d0))
    enddo
    s2 = s2 + f(b-dx/2.d0)
    s = (s1*2.d0 + s2*4.d0 + f(a) + f(b))*(dx/6.d0)
  end function simple_integrate


  function testf(x)
    COOP_REAL::x, testf
    testf = sin(x)
  end function testf
  
end module tmp


program searchzeq
  use coop_wrapper_utils
  use tmp
  !#define DOINT(f, a, b) simple_integrate(f, a, b)
#define DOINT(f, a, b) simple_integrate(f, a, b)  
  COOP_INT::i
  COOP_REAL::eps
  COOP_REAL::s1,s2,s
  print*, DOINT(testf, 0.d0, coop_pi)/2.d0-1.d0
  do i=1,100,3
     eps = sqrt(0.5d0*theta0*(i/100.d0))
     s1 = DOINT(intdt, 0.d0, theta0-2.d0*eps**2)  
     s2 = DOINT(intdx, 0.d0, eps)*4.d0
     s = (s1+s2)/coop_pi
     print*,i, s, 1.d0+theta0**2/16.d0+3.58d-3*theta0**4
  enddo
end program searchzeq

