program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT, parameter::n = 1024*12+193
  COOP_REAL:: p(n), q(n)
  type(coop_file)::fp
  COOP_INT::i
  do i=1, n
     p(i) = i
  enddo
  call fp%open("numbers.dat", 'u')
  read(fp%unit) q
  call fp%close()
  write(*,*) maxval(abs(p-q))
end program Test  
