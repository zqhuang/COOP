program test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_STRING str
  str = trim(coop_num2str(-0.722099998d8))
  write(*,*) str
  
end program test
