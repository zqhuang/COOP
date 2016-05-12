program test
  use coop_wrapper_firstorder
  use coop_ellipse_collapse_mod
  implicit none
#include "constants.h"
  type(coop_ellipse_collapse_params)::params
  logical::want_debug 
  COOP_INT, parameter::nz = 800
  COOP_REAL, parameter::zmax_check = 20.d0
  COOP_REAL, parameter::Omega_m = 0.3d0
  COOP_REAL,parameter::w = -1.d0
  COOP_REAL::F_pk, p_nu, e_nu, y(6), z(nz), a(nz), adynamic, z_collapse(3)
  COOP_INT::iz, iaxis
  call coop_get_command_line_argument(key = "debug", arg = want_debug, default = .false.)
  write(*,*) "========================================================"
  write(*,*) "Enter F_pk, e_nu, p_nu"
  call params%init(Omega_m = Omega_m, w = w)
  call coop_set_uniform(nz, z, zmax_check, 0.d0)
  a = 1.d0/(1.d0+z)
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
     do iz = 1, nz
        call params%evolve(adynamic, y, a(iz))
        if(want_debug)then
           write(*,"(4E16.7)") z(iz), y(1:3)/adynamic
        endif
        do iaxis  = 1, 3
           if(z_collapse(iaxis) .lt. 0.d0)then
              if(y(iaxis)/adynamic .le. params%collapse_a_ratio(iaxis)) & 
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
