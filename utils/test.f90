program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  !!this program find the zeros of J_n and the integrals between zeros
  type(coop_hankel_kernel)::hk
  COOP_REAL,parameter::a = 0.8d0  
  COOP_INT,parameter::n=50
  COOP_REAL::l(n), s(n)
  COOP_INT::i
  call coop_set_uniform(n, l, 0.1d0, 3.d0)
  call hk%init(0)
  do i = 1, n
     s(i) = hk%integrate(f, l(i))
     write(*,*) l(i), g(l(i)), s(i),  s(i) / g(l(i)) - 1.d0
  enddo

!!$  do i=2, hk%n-1
!!$     print*, hk%x(i),  hk%kernel(i)*2.d0/(hk%x(i+1)-hk%x(i-1)), besjn(0, hk%x(i))*hk%x(i)
!!$  enddo

contains

  function f(r)
    COOP_REAL::r, f
    f = exp(-0.5d0*a**2*r**2)
  end function f

  function g(k)
    COOP_REAL::k, g
    g = exp(-k**2/a**2/2.d0)/a**2
  end function g
  
end program Test  
