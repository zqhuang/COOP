program plot
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  call coop_file_decrypt("../../../doc/info/usr/usr.txt", "deusr.txt")
end program plot
