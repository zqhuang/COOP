program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_math_expression)::cme
  integer error
  COOP_REAL::res, ans
  res = -sin(-12.5 )*exp(-1.2-2.**3 / 5.2)-14./53.
  call coop_eval_math("$1 + $2", ans, (/ 14.d0, 5.2d0, 2.d0 /))
  print*, res, ans
end program Test
