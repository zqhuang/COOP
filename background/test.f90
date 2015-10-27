program test
  use coop_background_mod
  use coop_wrapper_utils
#include "constants.h"
  type(coop_function)::alpha
  COOP_INT,parameter::n=64
  COOP_REAL::lna(n)
  COOP_INT::i
  alpha = coop_de_alpha_constructor(0.1d0, "linear")
  call coop_set_uniform(n, lna, -1.d0, 0.d0)
  do i=1, n
     write(*,*) lna(i), alpha%eval(exp(lna(i)))
  enddo
end program test
