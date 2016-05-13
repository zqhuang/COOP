program test
  use coop_wrapper_firstorder
  use coop_ellipse_collapse_mod
  implicit none
#include "constants.h"
  type(coop_ellipse_collapse_params)::params
  logical::want_debug 
  COOP_INT, parameter::nz = 5000
  COOP_REAL, parameter::zmax_check = 10.d0
  COOP_REAL::Omega_m,w 
  COOP_REAL::F_pk, p_nu, e_nu, y(6), z(nz), a(nz), adynamic, z_collapse(3)
  COOP_INT::iz, iaxis
  call coop_get_command_line_argument(key = "debug", arg = want_debug, default = .false.)
  call coop_get_command_line_argument(key = "omm", arg = Omega_m, default = 0.3d0)
  call coop_get_command_line_argument(key = "w", arg = w, default = -1.d0)
  write(*,*) "========================================================"
  write(*,*) "./Col [-omm ...(0.3)] [-w ...(-1)] [-debug ...(F)]"
  write(*,*) "Example:"
  write(*,*) "./Col -omm 0.31 -w -0.9  (run model Omega_m = 0.31 , w = -0.9; debug off)"
  write(*,*) "./Col -omm 0.31 -debug T (run model Omega_m = 0.31, default w = -1; debug on)"
  write(*,*) "./Col   (run model default Omega_m = 0.3 , default w = -1; debug off)"
  write(*,*) "========================================================"
  write(*,*) "Enter F_pk, e_nu, p_nu"
  call params%init(Omega_m = Omega_m, w = w)
!  params%zinit = 500.d0   !!use this to test the robustness of zinit (result does not change as long as zinit >> 1)
!  params%collapse_a_ratio = (/ 0.02d0, 0.171d0, 0.171d0 /)  !!you can use this to check the 1.686 threshold 
  if(want_debug)then
     call coop_set_uniform(nz, a, 1.d0/(1.d0+ params%zinit-0.1d0), 1.d0)
  else
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
        write(*,"(7E16.7)") params%zinit, y(1:3)/adynamic, y(4:6)
     endif
     do iz = 1, nz
        call params%evolve(adynamic, y, a(iz))
        if(want_debug)then
           write(*,"(7E16.7)") z(iz), y(1:3)/adynamic, y(4:6)
        endif
        do iaxis  = 1, 3
           if(z_collapse(iaxis) .lt. 0.d0)then
              if(y(iaxis)/adynamic / params%collapse_a_ratio(iaxis) .lt. 1.001 .and. y(iaxis +3) .le. 0.d0) & 
                   z_collapse(iaxis) = z(iz)
           endif
        enddo
        if(all(z_collapse .ge. 0.d0))exit
     enddo
     do iaxis = 3, 1, -1
        if(z_collapse(iaxis) .ge. 0)then
           write(*,"(A, F12.2)")"Z_Collapse_"//COOP_STR_OF(iaxis)//" = ", z_collapse(iaxis)
        else
           write(*,"(A)")  "Axis_"//COOP_STR_OF(iaxis)//" is not collapsed."
        endif
     enddo
     write(*,*) "========================================================"
     write(*,*) "Enter F_pk, e_nu, p_nu (enter a negative F_pk to finish)"
  enddo
end program test
