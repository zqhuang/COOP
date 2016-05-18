program RunF
  use coop_wrapper_firstorder
  use coop_fisher_mod
  implicit none
#include "constants.h"
  type(coop_fisher)::fisher
  type(coop_file)::fp
  COOP_STRING::root, path
  COOP_INT::i, j
  if(iargc() .lt. 1) stop "Fisher input_file"
  call fisher%init(coop_InputArgs(1))
  root = trim(fisher%settings%value("root"))
  if(trim(root).eq."")then
     write(*,*) "you need to specify root in "//trim(coop_InputArgs(1))
     stop
  endif
  call fp%open(trim(root)//"_fisher.txt")  !!I do this first to test if the directory exists before spending a few minutes doing the fisher.
  call coop_prtsystime(.true.)
  call fisher%get_fisher()
  call coop_prtsystime()
  if(fisher%n_params_used .gt. 0)then
     write(fp%unit, "("//COOP_STR_OF(fisher%n_params_used)//"A16)") fisher%paramtable%key(fisher%ind_used)
     write(*,*) "saving fisher matrix to "//trim(root)//"_fisher.txt"
     do i=1, fisher%n_params_used
        write(fp%unit,"("//COOP_STR_OF(fisher%n_params_used)//"G16.7)") fisher%fisher(fisher%ind_used(i), fisher%ind_used) 
     enddo
     call fp%close()
     call fp%open(trim(root)//"_cov.txt")
     write(*,*) "saving covariance matrix to "//trim(root)//"_cov.txt"
     write(fp%unit, "("//COOP_STR_OF(fisher%n_params_used)//"A16)") fisher%paramtable%key(fisher%ind_used)
     write(fp%unit, "("//COOP_STR_OF(fisher%n_params_used)//"G16.7)") fisher%paramtable%val(fisher%ind_used)
     do i=1, fisher%n_params_used
        write(fp%unit,"("//COOP_STR_OF(fisher%n_params_used)//"G16.7)") fisher%cov(fisher%ind_used(i), fisher%ind_used) 
     enddo
     call fp%close()

     call fp%open(trim(root)//"_std.txt")
     write(*,*) "saving standard deviations to "//trim(root)//"_std.txt"
     do i=1, fisher%n_params_used
        write(fp%unit, "(A, G14.5, A5, G14.5)") trim(fisher%paramtable%key(fisher%ind_used(i)))//" = ", fisher%paramtable%val(fisher%ind_used(i)), " +/- ", sqrt(fisher%cov(fisher%ind_used(i), fisher%ind_used(i)))
        write(*, "(A, G14.5, A5, G14.5)") trim(fisher%paramtable%key(fisher%ind_used(i)))//" = ", fisher%paramtable%val(fisher%ind_used(i)), " +/- ", sqrt(fisher%cov(fisher%ind_used(i), fisher%ind_used(i)))
     enddo
     call fp%close()
     do i = 1, fisher%n_params_used
        do j = i+1, fisher%n_params_used
           if(fisher%cov(fisher%ind_used(i), fisher%ind_used(j))**2/fisher%cov(fisher%ind_used(i), fisher%ind_used(i))/fisher%cov(fisher%ind_used(j), fisher%ind_used(j)) .gt. 0.98)then
              write(*,*) "Warning: "//trim(fisher%paramtable%key(fisher%ind_used(i)))//" and  "//trim(fisher%paramtable%key(fisher%ind_used(j)))//" are strongly correlated (r="//COOP_STR_OF(fisher%cov(fisher%ind_used(i), fisher%ind_used(j))**2/fisher%cov(fisher%ind_used(i), fisher%ind_used(i))/fisher%cov(fisher%ind_used(j), fisher%ind_used(j)))//"). You may want to redefine the parameters to eliminate the degeneracy."
           endif
        enddo
     enddo
  else
     call fp%close()
     write(*,*) "Parameters are all fixed or unconstrained."
  endif

end program RunF
