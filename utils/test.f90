program plot
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  call coop_file_encrypt("../../../doc/info/usr/usr.txt", "usr.txt")
  call coop_file_decrypt("usr.txt", "deusr.txt")
end program plot
