program getbg
  use coop_wrapper_utils
  use coop_lattice_fields_mod
  implicit none
#include "constants.h"
  call coop_infbg%setup(nflds = 2, epsilon_end = 1.d0, f_ini = (/ 5.d0*coop_lattice_Mp, 6.d-6*coop_lattice_Mp /) )

  write(*,"(A)") "================================================================================="
  write(*, "(A, E16.7)") "Inflation lasted for :", -coop_infbg%lna(1), " efolds"
  write(*,"(A)") "================================================================================="
  write(*, "(A, E16.7)") "H/M_p at the end of inflation:", exp(coop_infbg%lnH(coop_infbg%nsteps))
  write(*,"(A)") "---------------------------------------------------------------------------------"    
  write(*, "(A)") "Field values (/M_p) at the end of inflation:"
  write(*, "("//COOP_STR_OF(coop_infbg%nflds)//"E16.7)") coop_infbg%f(coop_infbg%nsteps, :)/coop_lattice_Mp
  write(*,"(A)") "---------------------------------------------------------------------------------"      
  write(*, "(A)") "Time derivatives of field values (/M_p^2) at the end of inflation:"
  write(*, "("//COOP_STR_OF(coop_infbg%nflds)//"E16.7)") coop_infbg%fd(coop_infbg%nsteps, :)/coop_lattice_Mpsq
  write(*,"(A)") "---------------------------------------------------------------------------------"

end program Getbg
