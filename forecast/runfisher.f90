program RunF
  use coop_wrapper_firstorder
  use coop_fisher_mod
  implicit none
#include "constants.h"
  type(coop_fisher)::fisher
  type(coop_file)::fp
  COOP_STRING::root
  COOP_INT::i, j
  if(iargc() .lt. 1) stop "Fisher input_file"
  call fisher%init(coop_InputArgs(1))
  root = trim(fisher%settings%value("root"))
  if(trim(root).eq."")then
     write(*,*) "you need to specify root in "//trim(coop_InputArgs(1))
     stop
  endif

  call fisher%get_fisher()

  call fp%open(trim(root)//"_fisher.txt")
  write(fp%unit, "(A2, "//COOP_STR_OF(fisher%n_params_used)//"A16)") "# ", fisher%paramtable%key(fisher%ind_used)
  do i=1, fisher%n_params_used
     write(fp%unit,"("//COOP_STR_OF(fisher%n_params_used)//"G16.7)") fisher%fisher(fisher%ind_used(i), fisher%ind_used) 
  enddo
  call fp%close()


  call fp%open(trim(root)//"_cov.txt")
  write(fp%unit, "(A2, "//COOP_STR_OF(fisher%n_params_used)//"A16)") "# ", fisher%paramtable%key(fisher%ind_used)
  do i=1, fisher%n_params_used
     write(fp%unit,"("//COOP_STR_OF(fisher%n_params_used)//"G16.7)") fisher%cov(fisher%ind_used(i), fisher%ind_used) 
  enddo
  call fp%close()

  call fp%open(trim(root)//"_std.txt")
  do i=1, fisher%n_params_used
     write(fp%unit, "(A, G14.5, A5, G14.5)") trim(fisher%paramtable%key(fisher%ind_used(i)))//" = ", fisher%paramtable%val(fisher%ind_used(i)), " +/- ", sqrt(fisher%cov(fisher%ind_used(i), fisher%ind_used(i)))
     write(*, "(A, G14.5, A5, G14.5)") trim(fisher%paramtable%key(fisher%ind_used(i)))//" = ", fisher%paramtable%val(fisher%ind_used(i)), " +/- ", sqrt(fisher%cov(fisher%ind_used(i), fisher%ind_used(i)))
  enddo
  call fp%close()

end program RunF
