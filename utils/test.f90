program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  print*, coop_dir_exists("testdir/")
  
end program Test
