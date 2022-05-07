program plot
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  !!  call coop_file_decrypt("../../../doc/info/usr/usr.txt", "deusr.txt")
  if(coop_file_exists("deusr.txt"))then
     call coop_file_encrypt("deusr.txt", "usr.txt")
     !call system("scp usr.txt gw.cita.utoronto.ca:~/work/doc/info/usr/")  
     call system("mv usr.txt ../../../doc/info/usr/usr.txt")
     !call system("rm -f deusr.txt*")
  else
     print*, "Cannot find the file deusr.txt!"
  endif
end program plot
