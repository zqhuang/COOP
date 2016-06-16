program test
  use coop_wrapper_firstorder
  use coop_coupledDE_collapse_mod
  implicit none
#include "constants.h"
  type(coop_dictionary)::dict
  COOP_STRING::params_file
  type(coop_coupledDE_collapse_params)::params
  COOP_INT,parameter::na = 801
  COOP_REAL::zvir1
  COOP_REAL::F_pk, p_nu, e_nu
  COOP_INT::i
  COOP_REAL::a(na), x(3, na)
  type(coop_file)::fp
  COOP_STRING::output
  if(iargc().lt. 2)then
     write(*,*) "========================================================"
     write(*,*) "Syntax:"
     write(*,*) "./CDSolve params.ini"
     write(*,*) "========================================================"
  endif
  call coop_get_Input(1, params_file)
  call coop_load_dictionary(params_file, dict)
  call coop_dictionary_lookup(dict, "output", output)
  call coop_dictionary_lookup(dict, "collapse_ratio1", params%collapse_a_ratio(1), 0.178d0)
  call coop_dictionary_lookup(dict, "collapse_ratio2", params%collapse_a_ratio(2), 0.178d0)
  call coop_dictionary_lookup(dict, "collapse_ratio3", params%collapse_a_ratio(3), 0.178d0)
  call coop_dictionary_lookup(dict, "collapse_fpk", F_pk)
  call coop_dictionary_lookup(dict, "collapse_e", e_nu)
  call coop_dictionary_lookup(dict, "collapse_p", p_nu)

  params%collapse_a_ratio = max(params%collapse_a_ratio, 0.01d0)
  call params%init(dict, update_cosmology = .true.)
  call coop_set_uniform(na, a, 0.02d0, 1.d0)
  call params%get_solution(a, x)
  call fp%open(output)
  do i=1, na
     write(fp%unit, "(4E16.7)") a(i), x(:, i)
  enddo
  call fp%close()
  write(*,"(A)") "The solution is successfully written to "//trim(output)//"."
end program test
