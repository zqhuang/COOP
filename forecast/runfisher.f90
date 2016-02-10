program RunMC
  use coop_wrapper_firstorder
  use coop_fisher_mod
  implicit none
#include "constants.h"
  type(coop_fisher)::fisher
  type(coop_file)::fp
  COOP_STRING::fname
  COOP_INT::i
  if(iargc() .lt. 1) stop "Fisher input_file"
  call fisher%init(coop_InputArgs(1))
  fname = trim(fisher%settings%value("output"))
  if(trim(fname).eq."")then
     write(*,*) "you need to specify output in "//trim(coop_InputArgs(1))
     stop
  endif

  call fisher%get_fisher()

  call fp%open(fname)
  do i=1, fisher%n_params_used
     write(fp%unit,"("//COOP_STR_OF(fisher%n_params_used)//"G14.5)") fisher%cov(fisher%ind_used(i), fisher%ind_used) 
  enddo
  
  do i=1, fisher%n_params_used
     write(fp%unit, "(A, G14.5)") trim(fisher%paramtable%key(fisher%ind_used(i))), sqrt(fisher%cov(fisher%ind_used(i), fisher%ind_used(i)))
     write(*, "(A, G14.5)") trim(fisher%paramtable%key(fisher%ind_used(i))), sqrt(fisher%cov(fisher%ind_used(i), fisher%ind_used(i)))
  enddo
  call fp%close()

end program RunMC
