program bgtest
  use coop_wrapper_background
  implicit none
#include "constants.h"  
  type(coop_cosmology_background)::bg
  COOP_INT, parameter::n = 512
  COOP_REAL:: lna(n),a(n), Qcpl, tracking_n, phi, phidot, V, Vp, Vpp
  COOP_INT::i, index_CDM , index_DE
  type(coop_ode)::ode
  type(coop_file)::fp
  !!DE  parameters
  Qcpl = 0.1d0 !!coupling between DE and CDM
  tracking_n = 0.d0 !!  V \propto 1 / phi^n (n > 0)

  !!=============  set up background ===============  
  call bg%init(h=0.7d0)
  call bg%add_species(coop_baryon(0.046d0))
  call bg%add_species(coop_radiation(bg%Omega_radiation()))
  call bg%add_species(coop_neutrinos_massless(bg%Omega_massless_neutrinos_per_species()*(bg%Nnu())))
  call coop_background_add_coupled_DE(bg, Omega_c = 0.25d0, Q = Qcpl, tracking_n = tracking_n) !, dlnQdphi = 0.2d0, dUdphi = 0.d0)
  
  index_CDM = bg%index_of("CDM")
  index_DE = bg%index_of("Dark Energy")

  !!
  call bg%setup_background()
  
  !!=============== test energy conservation ===========
#define LNHUBBLE y(1)
#define DLNHUBBLE_DLNA yp(1)

  call coop_set_uniform(n, lna, -18.d0, 0.d0)
  a = exp(lna)
  call ode%init(1)
  call ode%set_initial_conditions(lna(1), (/ log(bg%Hratio(a(1))) /) )
  call fp%open("background_output.dat","w")
  write(fp%unit, "(20A19)") "ln a", "ln H", "ln sqrt(rho/3)",  "phi", "dot phi", "V",  "w_DE", "ln(rho_m/rho_m0)"
  do i=1, n
     if(i.gt.1)call ode%evolve(fcn, lna(i))
     phi = bg%species(index_DE)%DE_phi(a(i))
     phidot = bg%species(index_DE)%DE_phidot(a(i))
     V =  bg%species(index_DE)%DE_V(phi)
     write(fp%unit, "(20E19.9)") lna(i), ode%LNHUBBLE, log(bg%Hratio(a(i))), phi, phidot, V, bg%species(index_DE)%wofa(a(i)) , log(bg%species(index_CDM)%density_ratio(a(i))), (phidot**2/2.d0-V)/(phidot**2/2.d0+V)   !!the difference between the 2nd column and the 3rd column is the relative error in energy conservation
  enddo
  call fp%close()

  !!===============  use gnuplot  to plot the results ============
  !!test energy conservation:  
  !!    plot "background_output.dat" u 1:($2-$3) w l
  !!test asymptotic ln rho_m a^3
  !!    plot "background_output.dat" u 1:($8+3*$1) w l
  !!test w_DE
  !!    plot "background_output.dat" u 1:7 w l,  "background_output.dat" u 1:(($5**2/2-$6)/($5**2/2+$6)) 
contains

  
  subroutine fcn(n, lna, y, yp)
    COOP_INT n
    COOP_REAL lna, y(1), yp(1)
    DLNHUBBLE_DLNA = bg%HdotbyHsq(exp(lna))
  end subroutine fcn
  
end program bgtest
