program test
  use coop_wrapper_firstorder
  use coop_ellipse_collapse_mod
  implicit none
#include "constants.h"
  type(coop_ellipse_collapse_params)::params
  COOP_INT,parameter::na = 801
  COOP_REAL::Omega_m,w, Omega_k, h, zvir1, wa, epsilon_s
  COOP_REAL::F_pk, p_nu, e_nu
  COOP_INT::i
  COOP_REAL::a(na), x(8, na)
  type(coop_file)::fp
  COOP_STRING::output
  if(iargc().lt. 2)then
     write(*,*) "========================================================"
     write(*,*) "./Solve -fpk ... -e ... -p ... -out ... [-omm ...(0.3)] [-omk ...(0.)] [-h ...(0.7)] [-w ...(-1)] [-wa ...(0)] [-eps ...(0)] [-fr1 ...(0.178)] [-fr2 ... (0.178)] [-fr3 ... (0.178)] "
     write(*,*) "Examples:"
     write(*,*) "./Solve -fpk 0  -e 0 -p 0 -eps 0.25 -out testnull.dat"                         
     write(*,*) "./Solve -fpk 3  -e 0.1 -p 0.1 -omm 0.31 -out savex_lcdm.dat"                    
     write(*,*) "./Solve -fpk 2.5 -e 0.2 -p 0.1 -omm 0.31 -w -0.9 -out savex_w.dat"
     write(*,*) "./Solve -fpk 3  -e 0.1 -p 0.1 -omm 0.31 -w -0.9 -wa 0.1 -out savex_w0wa.dat"
     write(*,*) "./Solve -fpk 3  -e 0.1 -p 0.1 -omm 0.31 -eps 0.2 -out savex_eps.dat"
     write(*,*) "./Solve -fpk 1.69 -e 0 -p 0 -omm 0.29 -out savex_spherical.dat"
     write(*,*) "./Solve -fpk 3 -e 0.1 -p 0.02 -omk 0.02 -h 0.68  -out savex_open.txt "
     write(*,*) "========================================================"
     stop
  endif

  call coop_get_command_line_argument(key = "omm", arg = Omega_m, default = 0.3d0)
  call coop_get_command_line_argument(key = "omk", arg = Omega_k, default = 0.d0)
  call coop_get_command_line_argument(key = "h", arg = h, default = 0.7d0)
  call coop_get_command_line_argument(key = "w", arg = w, default = -1.d0)
  call coop_get_command_line_argument(key = "wa", arg = wa, default = 0.d0)
  call coop_get_command_line_argument(key = "eps", arg = epsilon_s, default = 0.d0)  
  call coop_get_command_line_argument(key = "fr1", arg = params%collapse_a_ratio(1), default = 0.178d0)
  call coop_get_command_line_argument(key = "fr2", arg = params%collapse_a_ratio(2), default = 0.178d0)
  call coop_get_command_line_argument(key = "fr3", arg = params%collapse_a_ratio(3), default = 0.178d0)
  call coop_get_command_line_argument(key = "fpk", arg = F_pk)
  call coop_get_command_line_argument(key = "e", arg = e_nu)
  call coop_get_command_line_argument(key = "p", arg = p_nu)
  call coop_get_command_line_argument(key = "out", arg = output)

  params%collapse_a_ratio = max(params%collapse_a_ratio, 2.d-3)

  call params%init(Omega_m = Omega_m, Omega_k = Omega_k, h = h, w = w, wa = wa, epsilon_s = epsilon_s, F_pk = F_pk, e_nu = e_nu, p_nu = p_nu)
  call coop_set_uniform(na, a, 0.03d0, 1.d0)
  call params%get_solution(a, x)
  call fp%open(output)
  write(fp%unit, "(9A16)")  "# a             ", " x1 ", " x2 ", " x3 ", " dot x1 ", " dot x2 ", " dot x3 ", " H a^{3/2} ", " D/a "
  do i=1, na
     write(fp%unit, "(9E16.7)") a(i), x(:, i)
  enddo
  call fp%close()
  write(*,"(A)") "The solution is successfully written to "//trim(output)//"."
end program test
