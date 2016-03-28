program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT::l, i
  COOP_REAL::x, amp, phase
  l = 100
  call coop_jl_setup_amp_phase(l)
  do  i = 80, 200
     x = i
     call coop_jl_get_amp_phase(l, x, amp, phase)
     print*, amp*cos(phase), coop_jl(l, x)
  enddo
end program Test  
