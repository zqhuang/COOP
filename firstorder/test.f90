program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  COOP_INT, parameter:: n = 128
  type(coop_cosmology_firstorder)::cosmology
  COOP_REAL lna(n), a(n), omegab(n), omegac(n), omegag(n), omeganu(n), omegal(n), omega_all
  COOP_INT::i
  type(coop_asy)::fig
  !!initialize cosmology
  call cosmology%set_standard_cosmology(Omega_b=0.05d0, Omega_c=0.265d0, h = 0.68d0, tau_re = 0.078d0, As = 2.195d-9, ns = 0.9655d0, nu_mass_ev = 0.06d0)
  call coop_set_uniform(n, lna, log(0.7d-5), 0.d0)
  a = exp(lna)
  do i=1, n
     omegab(i) = O0_BARYON(cosmology)%density(a(i))
     omegac(i) = O0_CDM(cosmology)%density(a(i))
     omegag(i) = O0_RADIATION(cosmology)%density(a(i))
     omeganu(i) = O0_NU(cosmology)%density(a(i)) + O0_MASSIVENU(cosmology)%density(a(i))
     omegal(i) = O0_DE(cosmology)%density(a(i))
     omega_all = omegab(i)  + omegac(i)  + omegag(i)  + omeganu(i) + omegal(i)
     omegab(i) = omegab(i)/omega_all
     omegac(i) = omegac(i)/omega_all
     omegag(i) = omegag(i)/omega_all
     omeganu(i) = omeganu(i)/omega_all
     omegal(i) = omegal(i)/omega_all
  enddo
  call fig%open("species.txt")
  call fig%init(xlabel = "$a$", ylabel = "energy fraction", xlog = .true., ylog = .true., ymin = 1.e-5, ymax = 1., xmin = real(a(1)), xmax = real(a(n)), doclip = .true.)
  call fig%curve(x = a, y = omegal, legend="Dark Energy", color = fig%color(1), linetype = fig%linetype(1), linewidth = 2.)  
  call fig%curve(x = a, y = omegac, legend="CDM", color = fig%color(2), linetype = fig%linetype(2), linewidth = 2.)
  call fig%curve(x = a, y = omegag, legend="photons", color = fig%color(3), linetype = fig%linetype(3), linewidth = 2.)
  call fig%curve(x = a, y = omeganu, legend="neutrinos", color = fig%color(4), linetype = fig%linetype(4), linewidth = 2.)  
  call fig%curve(x = a, y = omegab, legend="baryon", color = fig%color(5), linetype = fig%linetype(5), linewidth = 2.)


  call fig%legend(0.1, 0.45)
  call fig%close()
end program test
