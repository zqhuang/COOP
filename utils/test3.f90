module weird
#include "constants.h"
  use coop_wrapper_utils
  type dy
     real*8,dimension(:,:),allocatable::cov
     integer::dim
  end type dy

  type ddy
     type(dy),dimension(2)::m
  end type ddy
  
end module weird
program searchzeq
  use coop_wrapper_utils
  use weird
  integer i
  type(ddy)::myddy
!  allocate(myddy%m(2))
  read(*,*) myddy%m(1)%dim
  allocate(myddy%m(1)%cov(myddy%m(1)%dim, myddy%m(1)%dim))
  myddy%m(1)%cov = 0.
  do i=1, myddy%m(1)%dim
     myddy%m(1)%cov(i,i) = 1.+i
  enddo
  do i=1, myddy%m(1)%dim
     write(*,*) myddy%m(1)%cov(i,:)
  enddo
  print*
  call coop_matrix_inverse(myddy%m(1)%cov)
  do i=1, myddy%m(1)%dim
     write(*,*) myddy%m(1)%cov(i,:)
  enddo

  deallocate(myddy%m(1)%cov)
!  deallocate(myddy%m)
end program searchzeq

