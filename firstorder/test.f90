program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::fod
  type(coop_file)fp
  COOP_INT, parameter::nk=200
  COOP_INT ik
  COOP_REAL k(nk), Pk(nk), phi(nk), lnk(nk)
  call coop_set_uniform(nk, lnk, log(0.3d0), log(2.d3))
  k = exp(lnk)
  
  !!set cosmology
  call fod%Set_Planck_bestfit()
  !!print*, fod%zre
  !!if you want extended models
  !!call fod%set_standard_cosmology(Omega_b=0.047d0, Omega_c=0.26d0, h = 0.68d0, tau_re = 0.08d0, nu_mass_eV = 0.06d0, As = 2.15d-9, ns = 0.962d0, nrun = -0.01d0, r = 0.2d0, nt = -0.01d0, YHe = 0.25d0, Nnu = 3.d0)

  !!compute the scalar Cl's
  call fod%compute_source(0)

  print*, "sigma_8=",fod%sigma_R(r = 8.d0/fod%h()*fod%H0Mpc(), z=0.d0)

!!$


end program test
