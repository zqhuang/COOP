program Test
  use coop_wrapper_typedef
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT,parameter::n = 100
  COOP_REAL,dimension(:),allocatable::x
  type(coop_dynamic_array_integer)::inds
  call coop_MPI_init()
  allocate(x(0:n-1))
  x = 0.d0
  x(5:10) = 1.
  call inds%get_indices( x.gt.0.5d0 , start_index = 0)
  print*, inds%i
  call coop_MPI_Finalize()

end program Test
