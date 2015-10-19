program bgtest
  use coop_wrapper_background
  implicit none
#include "constants.h"
  COOP_REAL:: omegab = 0.049d0
  COOP_REAL:: omegac = 0.265d0
  type(coop_cosmology_background)::bg
  COOP_INT, parameter::nvars = 4
  COOP_INT, parameter::n = 512
  type(coop_function)::wp1, alpha_M
  COOP_INT::err
  type(coop_ode)::ode
  COOP_REAL::y(nvars), yp(nvars), a
  COOP_INT::index_de
  call wp1%init_polynomial( (/ 0.d0, 0.2d0 /) )
  call alpha_M%init_polynomial( (/ 0.d0, 0.d0, 0.d0, 0.d0 /) )
  call bg%set_alphaM(alpha_M)

  call bg%init(h=0.68d0)
  call bg%add_species(coop_baryon(omegab))
  call bg%add_species(coop_cdm(omegac))  
  call bg%add_species(coop_radiation(bg%Omega_radiation()))
  call bg%add_species(coop_neutrinos_massless(bg%Omega_massless_neutrinos_per_species()*(bg%Nnu())))

  call coop_background_add_EFT_DE_with_effective_w(bg, wp1, err)
  if(err .ne. 0) stop "cannot initialize EFT DE;"
  print*, bg%omega_k()
  
  call bg%setup_background()
  call ode%init(nvars)
  call ode%set_initial_conditions(0.d0, (/ bg%HdotbyHsq(1.d0), log(bg%Hratio(1.d0)), log(bg%species(5)%density(1.d0)), log(coop_Mpsq(1.d0)) /) )
  call ode%evolve(getderv, -1.d0)

  print*, ode%y(1),  ode%y(1)/bg%HdotbyHsq(exp(ode%x))-1.d0
  print*, ode%y(2), ode%y(2) - log(bg%Hratio(exp(ode%x)))
  print*, ode%y(3), ode%y(3) - log(bg%species(5)%density(exp(ode%x)))
  print*, ode%y(4), ode%y(4) - log(coop_Mpsq(exp(ode%x)))
contains

  subroutine getderv(n, lna, y, yp)
    COOP_INT::n
    COOP_REAL::lna, y(n), yp(n), a
    a = exp(lna)
    yp(1) = bg%HddbyH3(a) - bg%HdotbyHsq(a)**2 * 2.d0
    yp(2) = bg%HdotbyHsq(a)
    yp(3) = -3.d0*bg%species(5)%wp1effofa(a)
    yp(4) = coop_alphaM(a)
  end subroutine getderv

  
end program bgtest
