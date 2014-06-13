program derivative
  use coop_wrapper
  implicit none
#include "constants.h"
  COOP_INT, parameter::n  = 100000
  type(coop_function)::func_V
  COOP_REAL phi(n), V(n), x(n), y
  COOP_INT i
  type(coop_file) fp
  call fp%open("potential.txt")
  do i = 1, n
     read(fp%unit, *) phi(i), V(i)
  enddo
  phi =abs(phi)
  V = V + 9.d0
  if(any(phi .le. 0.d0) .or. (any(V.le.0.d0))) stop "phi, V<0"
  call func_V%init_nonUniform(x = phi, f = V, xlog = .true. , ylog = .true.)
  print*, func_V%xmin, func_V%xmax
  y = func_V%eval(1.d-4)
  print*, y
  y = func_V%derivative(1.d-4)
  print*, y
end program derivative
