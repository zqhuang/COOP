program test
  use coop_wrapper_firstorder
  use coop_ellipse_collapse_mod
  implicit none
#include "constants.h"
  type(coop_ellipse_collapse_params)::params
  COOP_REAL::Omega_m,w, Omega_k, h, zvir1
  COOP_REAL::F_pk, p_nu, e_nu
  call coop_get_command_line_argument(key = "omm", arg = Omega_m, default = 0.3d0)
  call coop_get_command_line_argument(key = "omk", arg = Omega_k, default = 0.d0)
  call coop_get_command_line_argument(key = "h", arg = h, default = 0.7d0)
  call coop_get_command_line_argument(key = "w", arg = w, default = -1.d0)
  call coop_get_command_line_argument(key = "fr1", arg = params%collapse_a_ratio(1), default = 0.18d0)
  call coop_get_command_line_argument(key = "fr2", arg = params%collapse_a_ratio(2), default = 0.18d0)
  call coop_get_command_line_argument(key = "fr3", arg = params%collapse_a_ratio(3), default = 0.18d0)
  params%collapse_a_ratio = max(params%collapse_a_ratio, 1.d-4)

  write(*,*) "========================================================"
  write(*,*) "./GetZ1 [-omm ...(0.3)] [-omk ...(0.)] [-h ...(0.7)] [-w ...(-1)] [-fr1 ...(0.18)] [-fr2 ... (0.18)] [-fr3 ... (0.18)] "
  write(*,*) "Examples:"
  write(*,*) "./GetZ1 -omm 0.31 -w -0.9"
  write(*,*) "./GetZ1 -omm 0.29"
  write(*,*) "./GetZ1 -omk 0.01 -h 0.68 -fr1 0.1 -fr2 0.1 -fr3 0.1"
  write(*,*) "========================================================"
  write(*,*) "Enter F_pk, e_nu, p_nu"
  call params%init(Omega_m = Omega_m, Omega_k = Omega_k, h = h, w = w)
  do
     read(*,*) F_pk, e_nu, p_nu
     if(F_pk .lt. 0.d0)then
        write(*,*) "negative F_pk: program done"
        exit
     endif
     call params%init(F_pk = F_pk, e_nu = e_nu, p_nu = p_nu)
     zvir1 = params%zvir1()
     write(*, "(A, F13.5)") "zvir1 = ", zvir1
     write(*,*) "========================================================"
     write(*,*) "Enter F_pk, e_nu, p_nu (enter a negative F_pk to finish)"
  enddo
end program test
