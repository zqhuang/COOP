program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  write(*,"(A)") trim(coop_string_strip_quotes(' "Tom"'))
end program Test
