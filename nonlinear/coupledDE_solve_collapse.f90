program test
  use coop_wrapper_firstorder
  use coop_coupledDE_collapse_mod
  implicit none
#include "constants.h"
#if DO_COUPLED_DE
  type(coop_dictionary)::dict
  COOP_STRING::params_file
  type(coop_coupledDE_collapse_params)::params
  COOP_INT,parameter::na = 801
  COOP_REAL::zvir1
  COOP_REAL::F_pk, p_nu, e_nu
  COOP_INT::i
  COOP_REAL::a(na), x(10, na)
  type(coop_file)::fp
  COOP_STRING::output

  if(iargc().lt. 1)then
     write(*,*) "========================================================"
     write(*,*) "Syntax:"
     write(*,*) "./CDSolve params.ini"
     write(*,*) "========================================================"
     stop
  endif
  call coop_get_Input(1, params_file)
  call coop_load_dictionary(params_file, dict)
  call coop_dictionary_lookup(dict, "output", output)
  call params%init(dict, update_cosmology = .true.)
  call coop_set_uniform(na, a, 0.02d0, 1.d0)
  call params%get_solution(a, x)
  call fp%open(output)
  write(fp%unit, "(11A16)")  "# a             ", " x1 ", " x2 ", " x3 ", " dot x1 ", " dot x2 ", " dot x3 ", " phi ", " dot phi ", " a^{-3/2}H ", " D/a "
  do i=1, na
     write(fp%unit, "(11E16.7)") a(i), x(:, i)
  enddo
  call fp%close()
  write(*,"(A)") "The solution is successfully written to "//trim(output)//"."
#endif
end program test
