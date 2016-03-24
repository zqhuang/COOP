program plot
  use coop_wrapper_utils
  implicit none
#include "constants.h"
!!  call coop_file_decrypt("../../../doc/info/usr/usr.txt", "deusr.txt")
  call coop_file_encrypt("deusr.txt", "usr.txt")
  call system("cp usr.txt ../../../doc/info/usr/usr.txt")
  call system("scp usr.txt gw.cita.utoronto.ca:~/work/doc/info/usr/")
  call system("rm -f deusr.txt")
end program plot
