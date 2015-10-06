program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_dictionary)::dict
  call coop_load_dictionary("tmp/test1.ini", dict)
  call dict%print()
end program Test
