program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  type(coop_function)::f
  COOP_INT:: i, j
  j = 2
  i = 1
  open(2, file='test.dat', access = 'stream')
  write(2) 1, 2, 3
  close(2)


end program Test  
