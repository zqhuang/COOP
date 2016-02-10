program RunMC
  use coop_wrapper_firstorder
  use coop_fisher_mod
  implicit none
#include "constants.h"
  type(coop_fisher)::fisher
  COOP_INT::i
  if(iargc() .lt. 1) stop "Fisher input_file"
  call fisher%init(coop_InputArgs(1))
  call fisher%get_fisher()
  do i=1, fisher%n_params_used
     write(*,"("//COOP_STR_OF(fisher%n_params_used)//"G14.5)") fisher%cov(fisher%ind_used(i), fisher%ind_used) 
  enddo
  
  do i=1, fisher%n_params_used
     write(*,*) trim(fisher%paramtable%key(fisher%ind_used(i))), sqrt(fisher%cov(fisher%ind_used(i), fisher%ind_used(i)))
  enddo

end program RunMC
