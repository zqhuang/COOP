program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT, parameter::n = 3
  COOP_REAL::invcov(n, n)
  COOP_INT::i, j
  invcov(1, :) = (/ 37.2026, 91.4793, -38.1707 /)

  invcov(2, :) = (/ 91.4793, 3556.39, -1803.65 /)

  invcov(3, :) = (/ -38.1707, -1803.65, 1445.78 /)

  call coop_sympos_inverse(n,n,invcov)
  write(*,"(A, F10.4)") "0.4412 +/- ", sqrt(invcov(3,3))


  invcov(1, :) = (/ 31.032, 77.773, -16.796 /)

  invcov(2, :) = (/  77.773, 2687.7, -1475.9 /)

  invcov(3, :) = (/ -16.796, -1475.9, 1323.0 /)

  call coop_sympos_inverse(n,n,invcov)
   write(*,"(A, F10.4)") "0.422 +/- ", sqrt(invcov(3,3))

end program Test  
