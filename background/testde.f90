program bgtest
  use coop_wrapper_background
  implicit none
#include "constants.h"  
  type(coop_cosmology_background)::bg
  COOP_INT, parameter::n = 512
  COOP_REAL:: lna(n),a(n), Qcpl, tracking_n, phi, phidot, V, Vp, Vpp, weff, Omegac, OmegaL, omegab, dlnQdphi, dUdphi, hub, rhom, rhode
  COOP_INT::i, index_CDM , index_DE
  type(coop_ode)::ode
  type(coop_file)::fp
  !!DE  parameters
  Qcpl = 0.3d0 !!coupling between DE and CDM
  tracking_n = 2.d0 !!  V \propto 1 / phi^n e^{U(phi)}   (n > 0)
  dlnQdphi = 0.1d0
  dUdphi = -1.3d0

  !!=============  set up background ===============
  omegab = 0.05d0
  omegac = 0.25d0

 
  call bg%init(h=0.7d0)
  call bg%add_species(coop_baryon(omegab))
  call bg%add_species(coop_radiation(bg%Omega_radiation()))
  call bg%add_species(coop_neutrinos_massless(bg%Omega_massless_neutrinos_per_species()*(bg%Nnu())))
  call coop_background_add_coupled_DE(bg, Omega_c = omegac, Q = Qcpl, tracking_n = tracking_n, dlnQdphi = dlnQdphi, dUdphi = dUdphi)
  
  index_CDM = bg%index_of("CDM")
  index_DE = bg%index_of("Dark Energy")

  !!
  call bg%setup_background()
  
  !!=============== test energy conservation ===========
#define LNHUBBLE y(1)
#define LNRHODE y(2)  
#define DLNHUBBLE_DLNA yp(1)
#define DLNRHODE_DLNA yp(2)
  
  call coop_set_uniform(n, lna, -17.d0, 0.d0)
  a = exp(lna)
  call ode%init(2)
  call ode%set_initial_conditions(lna(1), (/ log(bg%Hratio(a(1))),  log(bg%species(index_DE)%density(a(1))) /) )
  call fp%open("background_output.dat","w")
  write(fp%unit, "(20A19)") "#ln a", "ln H", "ln sqrt(rho/3)",  "phi", "dot phi", "V",  "w_DE", "ln(rho_m/rho_m0)"
  do i=1, n
     if(i.gt.1)call ode%evolve(fcn, lna(i))
     phi = bg%species(index_DE)%DE_phi(a(i))
     phidot = bg%species(index_DE)%DE_phidot(a(i))
     V =  bg%species(index_DE)%DE_V(phi)
     hub = exp(ode%LNHUBBLE)
     rhom = bg%species(index_CDM)%density(a(i))
     rhode =bg%species(index_DE)%density(a(i))
     if(a(i).gt. 0.2)then
        weff = - ( (bg%species(index_CDM)%dlnrhodlna(a(i))+3.d0)*rhom + (bg%species(index_DE)%dlnrhodlna(a(i))+3.d0)*rhode)/( Rhom + rhode - 3.d0*Omegac/a(i)**3 )/3.d0
     else
        weff = - ( (bg%species(index_CDM)%dlnrhodlna(a(i))+3.d0)*rhom + (bg%species(index_DE)%dlnrhodlna(a(i))+3.d0)*rhode)/rhode/3.d0
        weff = bg%species(index_DE)%wofa(a(i))
     endif
     write(fp%unit, "(20E19.9)") lna(i), a(i), ode%LNHUBBLE, log(bg%Hratio(a(i))), phi, V,  bg%species(index_DE)%wofa(a(i)) , log(bg%species(index_CDM)%density_ratio(a(i))), weff, rhom*a(i)**3/(3.d0*Omegac)
  enddo
  call fp%close()

  !!===============  use gnuplot  to plot the results ============
  !!test energy conservation:  
  !!    plot "background_output.dat" u 1:($3-$4) w l
  !!test asymptotic ln rho_m a^3
  !!    plot "background_output.dat" u 1:($8+3*$1) w l
  !!test w_DE, w_eff
  !!    plot "background_output.dat" u 2:7 w l,  "background_output.dat" u 2:9
contains

  
  subroutine fcn(n, lna, y, yp)
    COOP_INT n
    COOP_REAL lna, y(2), yp(2), a
    a = exp(lna)
    DLNHUBBLE_DLNA = bg%HdotbyHsq(a)
    DLNRHODE_DLNA = -3.d0*(1.d0+bg%species(index_DE)%wofa(a))
  end subroutine fcn
  
end program bgtest
