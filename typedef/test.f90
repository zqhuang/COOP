program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  integer, parameter::n = 3
  integer i
  do i=1,n
     call sleep(1)
     call cprint(".")
  enddo
  print*

end program Test
