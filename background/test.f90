program test
  use coop_background_mod
  use coop_wrapper_utils
#include "constants.h"
!!action = 1/2 \int (R + f(R)) ...
  type(coop_function)::fofR, Vofphi
  COOP_REAL::phi
  call fofR%init_powerlaw( c=(/ 0.1d0 /), alpha = (/ -2.d0 /) )
  call coop_convert_fofR_to_Vofphi(fofR, 2.1d0, Vofphi)
  do 
     read(*,*) phi
     if(phi .lt. 0.d0)exit
     print*, phi, Vofphi%eval(phi), Vofphi%derivative(phi), Vofphi%derivative2(phi)
  enddo
end program test
