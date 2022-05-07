program Test
  use coop_wrapper_utils
  use coop_lattice_fields_mod
  implicit none
#include "constants.h"
#include "lattice.h"
  type(coop_lattice_fields)::this
  type(coop_lattice_initial_power)::ini
  COOP_INT::i
  COOP_REAL::dt
  !!run ./BG to get the following initial background values 
  call coop_infbg%setup(nflds = 2, epsilon_end = 1.d0, f_ini = (/ 25.d0*coop_lattice_Mp, 0.d0 /) )
  call ini%alloc(nk = 20, nflds = 2)
  call ini%init(lnkmin = coop_infbg%lnHend - 5.d0, lnkmax = coop_infbg%lnHend+3.d0, is_diagonal = .true.)

  !!initialize   
  call this%init( n = 32, LH = 20.d0, phi = coop_infbg%fields(0.d0), pi = coop_infbg%dot_fields(0.d0), inipower= ini, use_conformal_time = .true.)
  !!set time step  
  dt = this%dx/10.d0
  !!set initial scale factor and Hubble parameter
  call this%set_pi_y()
  !!choose order for ODE solver (2, 4, or 6)
  this%ode_order = 6
  write(*,"(5A16)") " a ", "E_K/E_tot", " E_G/E_tot ", " Relative Error ", " a "
  write(*,"(5G16.7)") this%a, this%ke/(this%ke+this%ge+this%pe), this%ge/(this%ke+this%ge+this%pe), this%H**2*3.d0*coop_lattice_Mpsq / (this%ke + this%ge + this%pe) - 1.d0, this%scale_factor()
  !!evolve the fields
  do i = 1, 300
     call this%evolve(dt, 50)
     call this%set_energies()
     write(*,"(5G16.7)") this%a, this%ke/(this%ke+this%ge+this%pe), this%ge/(this%ke+this%ge+this%pe), this%H**2*3.d0*coop_lattice_Mpsq / (this%ke + this%ge + this%pe) - 1.d0, this%scale_factor()
  enddo
end program Test
