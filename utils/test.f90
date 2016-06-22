program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_STRING:: expr
  COOP_REAL::ans
  COOP_REAL,dimension(3)::vars = (/0.d0, 5.d0, 5.1d0 /)
  call coop_get_command_line_argument(key = "expr", arg = expr)
  print*, trim(expr)
!  expr = '($3-$2)/$2'
  call coop_eval_math(expr, ans, vars)
  print*, ans
end program Test  
