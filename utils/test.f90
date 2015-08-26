program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_dynamic_array_real)::raw, smooth
  COOP_STRING::input, output
  COOP_REAL::width
  COOP_INT::col
  COOP_INT::i
  if(iargc().lt.8)then
     write(*,*) "./SMOOTHDATA"
     stop
  endif
     
  call coop_get_command_line_argument(key = 'inp', arg = input)
  call coop_get_command_line_argument(key = 'out', arg = output)  
  call coop_get_command_line_argument(key = 'width', arg = width)
  call coop_get_command_line_argument(key = 'col', arg = col)    
  call raw%load_txt(trim(input))
  
  
end program Test
