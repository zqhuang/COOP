program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  COOP_REAL,dimension(:,:),allocatable::ss
  COOP_INT::l1, l2
  read(*,*) l1, l2
  allocate(ss(l1:l2, 2:10))
  write(*,*) lbound(ss, 1), ubound(ss,1), lbound(ss,2), ubound(ss,2)
end program test
