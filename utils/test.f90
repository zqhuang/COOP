program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_REAL,dimension(5)::vars = (/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0 /)
  COOP_REAL::ans
  call coop_eval_math('$1 + ($2 + $3)^2', ans, vars)
  print*, ans

end program Test  
