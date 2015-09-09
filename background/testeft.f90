program bgtest
  use coop_wrapper_background
  implicit none
#include "constants.h"  
  type(coop_cosmology_background)::bg
  COOP_INT, parameter::nvars = 1
  COOP_INT, parameter::n = 512
  type(coop_function)::wp1, alpha_M
  COOP_INT::err
  type(coop_ode)::ode
  COOP_REAL:: lna(n), lnH(n), omegab, omegac
  COOP_REAL::y(1), yp(1)
  COOP_INT::index_de
  omegab = 0.049d0
  omegac = 0.265d0
  call wp1%init_polynomial( (/ 0.d0, -0.2d0 /) )
  call alpha_M%init_polynomial( (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.2d0 /) )
  call bg%init(h=0.68d0)
  call bg%add_species(coop_baryon(omegab))
  call bg%add_species(coop_cdm(omegac))  
  call bg%add_species(coop_radiation(bg%Omega_radiation()))
  call bg%add_species(coop_neutrinos_massless(bg%Omega_massless_neutrinos_per_species()*(bg%Nnu())))
  call coop_background_add_EFT_DE(bg, wp1, alpha_M, err)  
  call bg%setup_background()
  call ode%init(nvars)
  call ode%set_initial_conditions(0.d0, (/ bg%HdotbyHsq(1.d0) /) )
  call ode%evolve(getderv, -1.d0)
  print*, ode%y(1), bg%HdotbyHsq(exp(ode%x)), ode%y(1)/bg%HdotbyHsq(exp(ode%x))-1.d0
  
contains

  subroutine getderv(n, lna, y, yp)
    COOP_INT::n
    COOP_REAL::lna, y(1), yp(1), a
    a = exp(lna)
    yp(1) = bg%HddbyH3(a) - bg%HdotbyHsq(a)**2 * 2.d0
  end subroutine getderv

  
end program bgtest
