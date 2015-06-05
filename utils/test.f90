program TestNpeak
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_SINGLE:: x
  call coop_get_commander_line_argument(key = 'xvalue', arg = x)
  print*, x
end program TestNpeak
