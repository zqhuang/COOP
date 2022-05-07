program LeastSQ
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT::nrows, col1, col2, ncols, i
  type(coop_dynamic_array_real)::arr
  COOP_REAL::k, b, r
  COOP_STRING::filename
  type(coop_file)::fp
  filename = coop_InputArgs(1)
  if(trim(filename) .eq. "")then
     write(*,*) "Enter number of data:"
     read(*,*) nrows
     ncols = 2
     col1 = 1
     col2 = 2
     call arr%init(nrows = nrows, ncols = ncols)
     do i=1, nrows
        write(*,*) "Enter data (x, y):"
        read(*,*) arr%f(i, :)
     enddo
  else
     call arr%load_txt(filename)
     if(arr%ncols .le. 1 .or. arr%nrows .le. 1)then
        write(*,*) "The file has "//COOP_STR_OF(arr%ncols)//" columns and "//COOP_STR_OF(arr%nrows)//" rows."
        write(*,*) "Cannot perform least square fit."
        stop
     endif
     call coop_get_command_line_argument(key = "x", arg = col1, default = 1)
     call coop_get_command_line_argument(key = "y", arg = col2, default = col1+1)
     if(col1 .gt. arr%ncols .or. col2 .gt. arr%ncols .or. col1.eq.col2)then
        write(*,*) "The file has "//COOP_STR_OF(arr%ncols)//" columns."
        write(*, *) "You want columns: "//COOP_STR_OF(col1)//" and "//COOP_STR_OF(col2)
        write(*,*) "Cannot perform least square fit."
        stop
     endif
  endif
  call coop_linear_least_square_fit(arr%nrows, arr%f(:, col1), arr%f(:, col2), k, b, r)
  write(*,"(A,G16.7)") "k = ", k
  write(*,"(A,G16.7)") "b = ", b
  write(*,"(A,G16.7)") "r = ", r
end program LeastSQ
