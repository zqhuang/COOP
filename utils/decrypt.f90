program plot
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  character yorn
  yorn = " "
  do while(yorn.ne."y" .and. yorn.ne."n" .and. yorn.ne."Y" .and. yorn.ne."N")
     print*, "Synchronize to the remote version? [Y|N]"
     read(*,*) yorn
  enddo
  if(yorn .eq. "Y" .or. yorn .eq. "y")then
     call system("scp gw.cita.utoronto.ca:~/work/doc/info/usr/usr.txt ./")
  else
     call system("cp ../../../doc/info/usr/usr.txt ./")
  endif
  call coop_file_decrypt("usr.txt", "deusr.txt")  
end program plot
