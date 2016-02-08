program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_STRING::str
  COOP_REAL::a
  !    print*, COOP_DATAPATH_I("%DATASETDIR%cmb/cmbdata_dataset%u.dat", 2)
  do 
     print*, "enter a number"
     read(*, "(A)") str
     if(coop_is_number(str))then
        read(str,*) a
        write(*,*) a
     else
        write(*,*) "not a number"
     endif
  enddo
end program Test
