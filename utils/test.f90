program TestNpeak
  use,intrinsic::iso_c_binding
  use coop_wrapper_utils
  implicit none
#include "constants.h"

  type(coop_asy_path)::path
  COOP_SINGLE::x0 = 1.2, y0 = 2.35
  COOP_SINGLE::p, a
  call path%append( x0 + 4., y0 + 3. )  
  call path%append( x0 + 4., y0 )  
  call path%append( x0, y0 )
  call path%close()
  call path%get_perimeter_and_area(1, p, a)
  write(*,*) p, a
end program TestNpeak
