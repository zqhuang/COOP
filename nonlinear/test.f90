program test
  use coop_wrapper_firstorder
  use coop_halofit_mod
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::fod
  type(coop_file)fp
  COOP_INT, parameter::nk=300
  COOP_INT ik
  COOP_REAL k(nk), Pl(nk), Pnl(nk), lnk(nk)

  call coop_set_uniform(nk, lnk, log(0.0001d0), log(2.5d0))
  k = exp(lnk)*(coop_SI_c/1.d5)
  !!set cosmology
  call fod%Set_Planck_bestfit()
  !!if you want extended models
  !!call fod%set_standard_cosmology(Omega_b=0.047d0, Omega_c=0.26d0, h = 0.68d0, tau_re = 0.08d0, nu_mass_eV = 0.06d0, As = 2.15d-9, ns = 0.962d0, nrun = -0.01d0, r = 0.2d0, nt = -0.01d0, YHe = 0.25d0, Nnu = 3.d0)

  !!compute the scalar Cl's
  call fod%compute_source(0)

  call coop_halofit_get_power(fod, 0.d0, nk, k, Pl, Pnl)
  k = k*(1.d5/coop_SI_c)
  Pl = Pl/k**3*(2.*coop_pi2)
  Pnl = Pnl/k**3*(2.*coop_pi2)

  call fp%open("pnl.txt", "w")
  do ik=1,nk
     write(fp%unit, "(20E16.7)") k(ik), Pl(ik), Pnl(ik)
  enddo
  call fp%close()


end program test
