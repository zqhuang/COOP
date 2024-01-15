program test
  use coop_wrapper_utils
  use coop_background_mod
  use coop_wrapper
  implicit none
#include "constants.h"
  logical,parameter::do_tracking = .true.
  type(coop_function)::V
  integer,parameter::n = 8192
  COOP_REAL, dimension(n)::phi_values, V_values
  logical,parameter::xlog = .true., ylog = .true.
  COOP_REAL,parameter::Omega_m = 0.3d0
  type(coop_arguments)::args
  type(coop_ode)::ode
  COOP_REAL::xend
  COOP_REAL::rhom_ini, phi_ini, V_ini, rhom, rhophi, epss, Gam, epsinf, phidot_ini, neff, meff, H_ini
  logical::search_eq
  type(coop_file)::fp
  COOP_REAL, parameter::Mpl = 1.d0
  COOP_INT i, j
  call coop_set_uniform(n, phi_values, 1.d-2, 5.d0, logscale = xlog)
  call args%init( r = (/ 1.d0 , 12.d0 /) )
  do i=1, n
     V_values(i) = potential(phi_values(i), args)
  enddo
  call V%init(n = n, xmin = phi_values(1), xmax = phi_values(n), f = V_values, xlog = xlog, ylog = ylog, check_boundary = .true.)
  rhom_ini = 1.d11*Mpl**2
  H_ini = sqrt(rhom_ini/(3.d0*Mpl**2))
  if(do_tracking)then
     !!random guess
     phi_ini = 0.05d0
     do j=1,10
        Gam = V%derivative2(phi_ini)*V%eval(phi_ini)/V%derivative(phi_ini)**2
        if(Gam .le. 1.d0) stop "bad initial guess"
        neff = 1.d0/(Gam - 1.d0)
        meff = 2.d0/(neff + 2.d0)
        phi_ini = (args%r(1)*neff/(meff*(meff+1.d0))/Mpl**2)**(meff/2.d0)*(2.d0/3.d0/H_ini)**meff
        phidot_ini = 1.5d0*meff*H_ini*phi_ini
        write(*, "(I5, 3E16.7)") j, phi_ini, V%eval(phi_ini)/rhom_ini, phidot_ini**2/2.d0/rhom_ini
     enddo
  else
     phi_ini = 0.01d0
     phidot_ini = 0.d0
  endif
  call ode%init(n = 2)
  call ode%set_initial_conditions( xini = 0.d0, yini = (/ phi_ini, phidot_ini /) )
  rhom = rhom_ini*exp(-3.d0*ode%x)
  rhophi = V%eval(ode%y(1))+ode%y(2)**2/2.d0
  search_eq = .true.
  epss = -1.d5
  call fp%open("w.txt", "w")
  write(fp%unit, "(10E16.7)") ode%x, (-V%eval(ode%y(1))+ode%y(2)**2/2.d0)/(V%eval(ode%y(1))+ode%y(2)**2/2.d0), rhophi/(rhom+rhophi)
  do while(rhom/rhophi .gt. Omega_m/(1.d0-Omega_m))
     xend = ode%x + 0.002d0
     call ode%evolve(eqs, xend)
     rhom = rhom_ini*exp(-3.d0*ode%x)
     rhophi = V%eval(ode%y(1))+ode%y(2)**2/2.d0
     if(search_eq)then
        if(rhom/rhophi .le. 0.501)then
           search_eq = .false.
           epss = (V%derivative(ode%y(1))/V%eval(ode%y(1)))**2 * Mpl**2/2.d0
           Gam = V%derivative2(ode%y(1))*V%eval(ode%y(1))/V%derivative(ode%y(1))**2
        endif
     endif
     write(fp%unit, "(10E16.7)") ode%x, (-V%eval(ode%y(1))+ode%y(2)**2/2.d0)/(V%eval(ode%y(1))+ode%y(2)**2/2.d0), rhophi/(rhom+rhophi)
  enddo
  call fp%close()
  print*, ode%y(1), V%eval(ode%y(1)), ode%y(2)**2/2.d0
  epsinf = 1.5d0/(1.d0+2.d0*(Gam-1.d0))
  print*, "epsilon_infinity = ", epsinf
  print*, "epsilon_s = ", epss


contains

  function potential(phi, args)
    COOP_REAL phi, potential
    type(coop_arguments)::args
    potential = args%r(1)/(phi/Mpl)**args%r(2) * exp((phi/Mpl)**2/2.d0)
!    potential = args%r(1) * (1.d0+cos(phi/(Mpl*args%r(2))))
  end function potential

  subroutine eqs(n, x, y, yp)
    COOP_INT n
    COOP_REAL x, y(n), yp(n), Hubble
    Hubble = sqrt((V%eval(y(1)) + y(2)**2/2.d0 + rhom_ini * exp(-3.d0*x))/(3.d0*Mpl**2))
    yp(1) = y(2)/Hubble
    yp(2) = -3.d0*y(2) - V%derivative(y(1))/Hubble
  end subroutine eqs


end program test
