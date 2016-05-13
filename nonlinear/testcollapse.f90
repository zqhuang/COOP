program test
  use coop_wrapper_firstorder
  use coop_ellipse_collapse_mod
  implicit none
#include "constants.h"
  type(coop_ellipse_collapse_params)::params
  logical::want_debug 
  COOP_INT::nz 
  COOP_REAL::zmax_check = 10.d0
  COOP_REAL::Omega_m,w, Omega_k, h
  COOP_REAL::F_pk, p_nu, e_nu, y(6), adynamic, z_collapse(3)
  COOP_INT::iz, iaxis
  COOP_REAL,dimension(:),allocatable::a, z!!theta is used for testing spherical solution
  call coop_get_command_line_argument(key = "debug", arg = want_debug, default = .false.)
  call coop_get_command_line_argument(key = "omm", arg = Omega_m, default = 0.3d0)
  call coop_get_command_line_argument(key = "omk", arg = Omega_k, default = 0.d0)
  call coop_get_command_line_argument(key = "h", arg = h, default = 0.7d0)
  call coop_get_command_line_argument(key = "w", arg = w, default = -1.d0)
  call coop_get_command_line_argument(key = "zini", arg = params%zinit, default = 200.d0)
  if(params%zinit .lt. 19.99d0)stop "zini must be >20"
  call coop_get_command_line_argument(key = "fr1", arg = params%collapse_a_ratio(1), default = 0.18d0)
  call coop_get_command_line_argument(key = "fr2", arg = params%collapse_a_ratio(2), default = 0.18d0)
  call coop_get_command_line_argument(key = "fr3", arg = params%collapse_a_ratio(3), default = 0.18d0)
  params%collapse_a_ratio = max(params%collapse_a_ratio, 1.d-4)

  write(*,*) "========================================================"
  write(*,*) "./Col [-omm ...(0.3)] [-omk ...(0.)] [-h ...(0.7)] [-w ...(-1)] [-fr1 ...(0.18)] [-fr2 ... (0.18)] [-fr3 ... (0.18)] [-zini ... (100)] [-debug ...(F)] "
  write(*,*) "Examples:"
  write(*,*) "./Col -omm 0.31 -w -0.9  (Omega_m = 0.31 , Omega_k = 0, w = -0.9; debug off)"
  write(*,*) "./Col -omm 0.31 -debug T (Omega_m = 0.31,  Omega_k = 0, w = -1; debug on)"
  write(*,*) "./Col -omk 0.01  (Omega_m = 0.3 , Omega_k=0.01, w = -1; debug off)"
  write(*,*) "========================================================"
  write(*,*) "Enter F_pk, e_nu, p_nu"
  call params%init(Omega_m = Omega_m, Omega_k = Omega_k, h = h, w = w)
  if(want_debug)then
     nz = 200     
     allocate(a(nz), z(nz))
     call coop_set_uniform(nz, a, 1.d0/(1.d0+ params%zinit), 1.d0)
  else
     nz = 10000     
     zmax_check = 10.d0
     allocate(a(nz), z(nz))
     call coop_set_uniform(nz, a, 1.d0/(1.d0+min(zmax_check, params%zinit-0.1d0)), 1.d0)
  endif
  z = 1.d0/a - 1.d0

  do
     z_collapse = -1.d0
     read(*,*) F_pk, e_nu, p_nu
     if(F_pk .lt. 0.d0)then
        write(*,*) "negative F_pk: program done"
        exit
     endif
     call params%init(F_pk = F_pk, e_nu = e_nu, p_nu = p_nu )
     call params%set_initial_conditions(y)
     adynamic = 1.d0/(1.d0+params%zinit)
     if(want_debug)then
        write(*,"(7E16.7)") 1.d0/(1.d0+params%zinit), y(1:3)/adynamic, y(1:3)
     endif
     do iz = 1, nz
        call params%evolve(adynamic, y, a(iz))
        if(want_debug)then
           write(*,"(7E16.7)") a(iz), y(1:3)/adynamic, y(1:3)
        endif
        do iaxis  = 1, 3
           if(z_collapse(iaxis) .lt. 0.d0)then
              if(y(iaxis)/adynamic / params%collapse_a_ratio(iaxis) .lt. 1.001 .and. y(iaxis +3) .lt. 0.d0) & 
                   z_collapse(iaxis) = z(iz)
           endif
        enddo
        if(all(z_collapse .ge. 0.d0))exit
     enddo
     do iaxis = 3, 1, -1
        if(z_collapse(iaxis) .ge. 0)then
           write(*,"(A, F12.3)")"Z_Collapse_"//COOP_STR_OF(iaxis)//" = ", z_collapse(iaxis)
        else
           write(*,"(A)")  "Axis_"//COOP_STR_OF(iaxis)//" is not collapsed."
        endif
     enddo
     write(*,*) "========================================================"
     write(*,*) "Enter F_pk, e_nu, p_nu (enter a negative F_pk to finish)"
  enddo
end program test
