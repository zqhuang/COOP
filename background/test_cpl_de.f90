program bgtest
  use coop_wrapper_background
  implicit none
#include "constants.h"
  COOP_REAL, parameter::w0 = 0.1d0, wa = -0.1d0
  COOP_REAL, parameter::Q0 = 0.1d0, Qa = -0.1d0
  COOP_REAL, parameter::Omegac = 0.25d0, Omegab = 0.05d0
  COOP_REAL::lna
  type(coop_cosmology_background)::bg
  COOP_INT::i, index_CDM , index_DE, err
  type(coop_ode)::ode
  type(coop_function)::fwp1, fQ, Vofphi, intQofphi
  COOP_REAL::norm
#if DO_COUPLED_DE  
!  call fwp1%init_polynomial( (/ 0.1d0, 0.05d0 /) )
!  call fQ%init_polynomial( (/ 0.1d0, 0.05d0 /) )
  call Vofphi%init_powerlaw( c = (/ 2.d0, 0.1d0 /), alpha = (/ 0.d0, -0.1d0 /), name = "V(phi)")

  call intQofphi%init_polynomial( (/ 0.d0, 0.d0/) )

  call bg%init(h=0.7d0)
  call bg%add_species(coop_baryon(omegab))
  call bg%add_species(coop_radiation(bg%Omega_radiation()))
  call bg%add_species(coop_neutrinos_massless(bg%Omega_massless_neutrinos_per_species()*(bg%Nnu())))
!  call coop_background_add_coupled_DE(bg, Omega_c = omegac, fwp1 = fwp1, fQ = fQ, err = err)
  call coop_background_add_coupled_DE_with_potential(bg, Omega_c = omegac, Vofphi = Vofphi, intQofphi = intQofphi, norm = norm, err = err)

  if(err .ne. 0)then
     print*, err
  else
     index_CDM = bg%index_of("CDM")
     index_DE = bg%index_of("Dark Energy")
     call ode%init(n = 3)
     lna = -5.d0
     call ode%set_initial_conditions(xini = lna, yini = (/ log(bg%Hratio(exp(lna))), log(bg%species(index_DE)%density(exp(lna))),  log(bg%species(index_CDM)%density(exp(lna))) /) )
     call ode%evolve(lnHeq, xend = 0.d0)
     print*, ode%y(1) - log(bg%Hratio(exp(ode%x)))
     print*, ode%y(2) - log(bg%species(index_DE)%density(exp(ode%x)))
     print*, ode%y(3) - log(bg%species(index_CDM)%density(exp(ode%x)))
     call coop_asy_plot_function(bg%species(index_DE)%fwp1, "wp1.txt")
     call coop_asy_plot_function(bg%species(index_DE)%fwp1eff, "wp1eff.txt")     
  endif
#else
  write(*,*) "Coupled DE disabled"
  stop "You need to enable it in configure.in (set DARK_ENERGY_MODEL = COUPLED_DE)"
  
#endif  
  
contains

  subroutine lnHeq(n, lna, y, yp)
    COOP_INT::n
    COOP_REAL y(n), yp(n), lna, a
    a = exp(lna)
    yp(1) = bg%HdotbyHsq(a)
    yp(2) = -3.d0*bg%species(index_DE)%wp1effofa(a)
    yp(3) = -3.d0*bg%species(index_CDM)%wp1effofa(a)
  end subroutine lnHeq
  
end program bgtest
