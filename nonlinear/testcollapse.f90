program test
  use coop_wrapper_firstorder
  use coop_ellipse_collapse_mod
  implicit none
#include "constants.h"
  type(coop_ellipse_collapse_params)::params
  COOP_INT::nz 
  COOP_REAL::zmax_check = 10.d0
  COOP_REAL::Omega_m,w, Omega_k, h, zvir1
  COOP_REAL::F_pk, p_nu, e_nu, a
  COOP_INT::i
  call coop_get_command_line_argument(key = "omm", arg = Omega_m, default = 0.3d0)
  call coop_get_command_line_argument(key = "omk", arg = Omega_k, default = 0.d0)
  call coop_get_command_line_argument(key = "h", arg = h, default = 0.7d0)
  call coop_get_command_line_argument(key = "w", arg = w, default = -1.d0)
  call coop_get_command_line_argument(key = "fr1", arg = params%collapse_a_ratio(1), default = 0.18d0)
  call coop_get_command_line_argument(key = "fr2", arg = params%collapse_a_ratio(2), default = 0.18d0)
  call coop_get_command_line_argument(key = "fr3", arg = params%collapse_a_ratio(3), default = 0.18d0)
  call coop_get_command_line_argument(key = "acc", arg = params%accuracy, default = 1.d-2)
  params%collapse_a_ratio = max(params%collapse_a_ratio, 1.d-4)

  write(*,*) "========================================================"
  write(*,*) "./Col [-omm ...(0.3)] [-omk ...(0.)] [-h ...(0.7)] [-w ...(-1)] [-fr1 ...(0.18)] [-fr2 ... (0.18)] [-fr3 ... (0.18)] [-acc ...(0.01)]"
  write(*,*) "Examples:"
  write(*,*) "./Col -omm 0.31 -w -0.9  (Omega_m = 0.31 , Omega_k = 0, h=0.7, w = -0.9;)"
  write(*,*) "./Col -omm 0.29  (Omega_m = 0.29,  Omega_k = 0, h=0.7, w = -1; )"
  write(*,*) "./Col -omk 0.01 -h 0.68 -acc 0.001 (Omega_m = 0.3 , Omega_k=0.01, h=0.68, w = -1; accuracy 0.001)"
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
     write(*,*) "lambda = ", params%lambda
     zvir1 = params%zvir1()
     write(*, "(A, F12.2)") "zvir1 = ", zvir1
     write(*,*) "========================================================"
     write(*,*) "Enter F_pk, e_nu, p_nu (enter a negative F_pk to finish)"
  enddo
end program test
