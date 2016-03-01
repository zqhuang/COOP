program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT::i
  COOP_REAL::vars(2), ans
  vars(1) = 0.5d0
  vars(2) = coop_pi
  call coop_eval_math("Gamma($1)/sqrt($2)+log10(100.)",  ans, vars)
  print*, ans
contains

  function gaussian(x)
    COOP_REAL::x, gaussian
    gaussian = exp(-x**2/2.d0)
  end function gaussian
end program Test
