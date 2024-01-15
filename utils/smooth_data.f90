program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_dynamic_array_real)::raw, smooth
  COOP_STRING::input, output
  COOP_REAL::width, sumw, w
  COOP_INT::col
  COOP_INT::i, j
  if(iargc().lt.8)then
     write(*,*) "Syntax:"
     write(*,*) "./SMOOTH -inp DATAFILE -out OUTPUT -width WIDTH -col COLUMN"
     stop
  endif
  call coop_get_command_line_argument(key = 'inp', arg = input)
  call coop_get_command_line_argument(key = 'out', arg = output)  
  call coop_get_command_line_argument(key = 'width', arg = width)
  call coop_get_command_line_argument(key = 'col', arg = col)    
  call raw%load_txt(trim(input))
  if(col .gt. raw%ncols) stop "data file does not contain enough columns"  
  smooth = raw
  do i=1, raw%nrows
     sumw = 0.d0
     smooth%f(i, col) = 0.d0     
     do j = 1, raw%nrows
        w = exp(-((i-j)/width)**2)
        sumw = sumw + w
        smooth%f(i, col) = smooth%f(i, col) + raw%f(j, col)*w
     enddo
     smooth%f(i, col) = smooth%f(i, col)/sumw
  enddo
  call smooth%dump_txt(trim(output))
  
end program Test
