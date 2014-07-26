program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  integer i
  type(coop_dictionary)::dict
  call coop_load_dictionary("test.txt", dict, col_key = 2, col_value = 3)
  call dict%print
end program Test
