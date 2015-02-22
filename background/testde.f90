program bgtest
  use coop_wrapper_background
  implicit none
#include "constants.h"  
  type(coop_cosmology_background)::bg
  COOP_INT, parameter::n = 256
  COOP_REAL:: a(n), lna(n), ln_rho_cdm(n), phi(n), lnH(n), phidot(n), V(n), KE(n)
  COOP_INT::i, index_CDM , index_DE
  type(coop_file)::fp
  call bg%init(h=0.68d0)
  call bg%add_species(coop_baryon(0.046d0))
  call bg%add_species(coop_radiation(bg%Omega_radiation()))
  call bg%add_species(coop_neutrinos_massless(bg%Omega_massless_neutrinos_per_species()*(bg%Nnu())))
  call coop_background_add_coupled_DE(bg, Omega_c = 0.25d0, Q = 0.4d0, tracking_n = 1.d0)
  index_CDM = bg%index_of("CDM")
  index_DE = bg%index_of("Dark Energy")
  call bg%setup_background()
  call coop_set_uniform(n, lna, -10.d0, 0.d0)
  a = exp(lna)
  call fp%open("background_output.dat","w")
  do i=1, n
     ln_rho_cdm(i) = log(bg%species(index_CDM)%density_ratio(a(i)))
     phi(i) = bg%species(index_DE)%DE_phi(a(i))
     phidot(i) = bg%species(index_DE)%DE_phidot(a(i))
     V(i) = bg%species(index_DE)%DE_V(phi(i))
     KE(i) = phidot(i)**2/2.d0
     write(fp%unit, "(10G15.6)") lna(i), ln_rho_cdm(i)+3.d0*lna(i), phi(i), phidot(i), V(i), bg%species(index_DE)%wofa(a(i))
  enddo
  call fp%close()
end program bgtest
