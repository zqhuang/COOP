program TestNpeak
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_REAL x
  call coop_MPI_init()
  open(10, file="test.dat", form="unformatted")
  write(10) 1.d0, 2.d0, 3.d0
  close(10)
  open(10, file="test.dat", form="unformatted", access="append")
  write(10) 4.d0, 5.d0
  close(10)
  open(10, file="test.dat", form="unformatted")
  read(10) x
  write(*,*) x
  close(10)
  call coop_MPI_finalize()
  
end program TestNpeak
