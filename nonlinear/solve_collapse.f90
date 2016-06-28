program test
  use coop_wrapper_firstorder
  use coop_ellipse_collapse_mod
  implicit none
#include "constants.h"
  type(coop_ellipse_collapse_params)::params
  COOP_INT,parameter::na = 801
  COOP_REAL::Omega_m,w, Omega_k, h, zvir1
  COOP_REAL::F_pk, p_nu, e_nu
  COOP_INT::i
  COOP_REAL::a(na), x(6, na)
  type(coop_file)::fp
  COOP_STRING::output
  if(iargc().lt. 2)then
     write(*,*) "========================================================"
     write(*,*) "./Solve -fpk ... -e ... -p ... -out ... [-omm ...(0.3)] [-omk ...(0.)] [-h ...(0.7)] [-w ...(-1)] [-fr1 ...(0.18)] [-fr2 ... (0.18)] [-fr3 ... (0.18)] "
     write(*,*) "Examples:"
     write(*,*) "./Solve -fpk 2.5 -e 0.2 -p 0.1 -omm 0.31 -w -0.9 -out savex.dat"
     write(*,*) "./Solve -fpk 1.6865 -e 0 -p 0 -omm 0.29 -out solution.dat"
     write(*,*) "./Solve -fpk 3 -e 0.1 -p 0.02 -omk 0.01 -h 0.68  -out mysol.txt "
     write(*,*) "========================================================"
     stop
  endif

  call coop_get_command_line_argument(key = "omm", arg = Omega_m, default = 0.3d0)
  call coop_get_command_line_argument(key = "omk", arg = Omega_k, default = 0.d0)
  call coop_get_command_line_argument(key = "h", arg = h, default = 0.7d0)
  call coop_get_command_line_argument(key = "w", arg = w, default = -1.d0)
  call coop_get_command_line_argument(key = "fr1", arg = params%collapse_a_ratio(1), default = 0.18d0)
  call coop_get_command_line_argument(key = "fr2", arg = params%collapse_a_ratio(2), default = 0.18d0)
  call coop_get_command_line_argument(key = "fr3", arg = params%collapse_a_ratio(3), default = 0.18d0)
  call coop_get_command_line_argument(key = "fpk", arg = F_pk)
  call coop_get_command_line_argument(key = "e", arg = e_nu)
  call coop_get_command_line_argument(key = "p", arg = p_nu)
  call coop_get_command_line_argument(key = "out", arg = output)

  params%collapse_a_ratio = max(params%collapse_a_ratio, 2.d-3)

  call params%init(Omega_m = Omega_m, Omega_k = Omega_k, h = h, w = w,F_pk = F_pk, e_nu = e_nu, p_nu = p_nu)
  call coop_set_uniform(na, a, 0.03d0, 1.d0)
  call params%get_solution(a, x)
  call fp%open(output)
  write(fp%unit, "(7A16)")  "# a             ", " x1 ", " x2 ", " x3 ", " dot x1 ", " dot x2 ", " dot x3 "
  do i=1, na
     write(fp%unit, "(7E16.7)") a(i), x(:, i)
  enddo
  call fp%close()
  write(*,"(A)") "The solution is successfully written to "//trim(output)//"."
end program test
