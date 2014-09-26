program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_function)fc
  call coop_file_load_function("deltaN-LUT-1.875", 1, 2, fc,.false.)
  print*, fc%eval(0.197014d-5)
end program Test
