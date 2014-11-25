program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_REAL n, chi2
  COOP_STRING::inline
  if(iargc() .ge. 2)then
     inline = coop_InputArgs(1)
     read(inline,*) n
     inline = coop_InputArgs(2)
     read(inline,*) chi2
  else
     write(*,*) "Enter n, chi^2 per dof "
     read(*,*) n, chi2
  endif
  print*, coop_IncompleteGamma(n/2.d0, chi2*n/2.d0)/gamma(n/2.d0)
end program Test
