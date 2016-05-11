program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  type(coop_2dfunction)::f2d
  COOP_REAL::x, y
  call f2d%init(f = f, nx = 100, ny = 100, xmin = -4.d0, xmax = 4.d0, ymin = -4.d0, ymax = 4.d0, zlog = .true.)
  do
     read(*,*) x, y
     if(x .lt. -1.e30 .or. x .gt. 1.e30 .or. y .le. -1.e30 .or. y.ge. 1.e30)exit
     print*, f2d%eval(x, y), f(x, y)
  enddo
contains

  function f(x, y)
    COOP_REAL::f, x, y
    f = exp(-x**2-y**2)*(x**2+1.d0)*(y**4+1.d0)
  end function f

end program Test  
