program Test
  use coop_wrapper_utils
  use coop_lattice_fields_mod
  implicit none
#include "constants.h"
#include "lattice.h"
  type(coop_lattice_fields)::this
  COOP_INT::i
  COOP_REAL::dt, phi(2), pi(2), phi_sigma2(2)
  phi  =  (/ coop_lattice_Mp, 0.d0 /)
  pi = (/ -coop_lattice_Mpsq*1.d-8, 0.d0 /)
  phi_sigma2 = (/ coop_lattice_Mpsq*1.d-14,  coop_lattice_Mpsq*1.d-14 /)
  call this%init( n = 32, LH = 10.d0, phi = phi, pi = pi, mu = 1.d0, phi_sigma2 = phi_sigma2)
  dt = this%dx/20.d0
  call this%set_pi_y()
  this%ode_order = 6.d0
  print*, this%a, this%ge/(this%ke+this%ge+this%pe), this%H**2*3.d0*coop_lattice_Mpsq / (this%ke + this%ge + this%pe) - 1.d0
  do i = 1, 100
     call this%evolve(dt, 50)
     call this%set_pi_y()
     print*, this%a, this%ge/(this%ke+this%ge+this%pe), this%H**2*3.d0*coop_lattice_Mpsq / (this%ke + this%ge + this%pe) - 1.d0
  enddo
end program Test
