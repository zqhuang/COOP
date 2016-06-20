program bgtest
  use coop_wrapper_background
  implicit none
#include "constants.h"
  COOP_REAL, parameter::Omega_c = 0.2442857d0, Omega_b = 0.045347
  COOP_REAL,parameter::Omega_m = omega_c+omega_b
  COOP_REAL,parameter::Omega_Lambda = 1.d0 - Omega_m
  COOP_REAL::lna  
  type(coop_cosmology_background)::bg
  COOP_INT::i, index_CDM , index_DE, err
  type(coop_ode)::ode
  type(coop_function)::fwp1, fQ, Vofphi, intQofphi
  type(coop_function)::fofR
  COOP_REAL::Lambda,  n, c2

#if DO_COUPLED_DE  

!  call fwp1%init_polynomial( (/ 0.1d0, 0.05d0 /) )
!  call fQ%init_polynomial( (/ 0.1d0, 0.05d0 /) )
!  call Vofphi%init_powerlaw( c = (/ 2.d0, 0.1d0 /), alpha = (/ 0.d0, -0.1d0 /), name = "V(phi)")
  !! action = 1/2 \int [R - 2 Lambda + f(R) ]...
  !!f(R) =   2 Lambda /(c2 * R^n + 1)
  Lambda = 2.1d0
  n = 1.d0
  c2 = n/Lambda**n/0.15

  call fofR%init_rational( c_up =(/ 2*Lambda /), alpha_up = (/ 0.d0 /), c_down = (/ c2, 1.d0 /),  alpha_down = (/ n, 0.d0 /) )
  call coop_convert_fofR_to_Vofphi(fofR, Lambda, Vofphi)

!  call coop_asy_plot_function(Vofphi, "Vofphi.txt")

  call intQofphi%init_polynomial( (/ 0.d0, 0.1d0 /) )

  call bg%init(h=0.7d0)
  call bg%add_species(coop_radiation(bg%Omega_radiation()))
  call bg%add_species(coop_neutrinos_massless(bg%Omega_massless_neutrinos_per_species()*(bg%Nnu())))
!  call coop_background_add_coupled_DE(bg, Omega_c = omega_c, fwp1 = fwp1, fQ = fQ, err = err)
  call coop_background_add_coupled_DE_with_potential(bg, Omega_c = omega_c, omega_b = omega_b, Vofphi = Vofphi, intQofphi = intQofphi,  err = err)


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
     
     print*, O0_DE(bg)%cplde_phi_lna%eval(0.d0)
     print*, O0_BARYON(bg)%rhoa3_ratio(1.d-10),   O0_CDM(bg)%rhoa3_ratio(1.d-10)
     print*, O0_BARYON(bg)%rhoa3_ratio(1.d-1),   O0_CDM(bg)%rhoa3_ratio(1.d-1)
     print*, O0_BARYON(bg)%rhoa3_ratio(0.5d0),   O0_CDM(bg)%rhoa3_ratio(0.5d0)

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
