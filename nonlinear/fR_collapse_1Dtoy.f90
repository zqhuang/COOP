program fR1d
  use coop_wrapper_firstorder
  use fR1d_mod
  implicit none
#include "constants.h"
  type(coop_fr1d_obj)::halo, haloini
  type(coop_asy)::fig, figphi
  COOP_INT::i, nsteps
  COOP_REAL::delta_ini, z_ini, delta_0, Omega_m, dlna, r_halo
  COOP_REAL,parameter::step = 5.d-8, fr = 1.d4
  nsteps = 50000
  Omega_m = 0.3d0
  z_ini = 99.d0
  delta_0 = 1.8d0
  r_halo = 0.001d0
  write(*, "(A, E15.4)") "M = ",  0.5d0 * r_halo**3 * coop_SI_PlanckMass * coop_SI_hbyH0/ coop_SI_PlanckTime / coop_SI_Msun / 0.7d0
  delta_ini = delta_0 * coop_Growth_fitting(Omega_m, -1.d0, z_ini) / coop_Growth_fitting(Omega_m, -1.d0, 0.d0)
  call halo%init(Omega_m = Omega_m, nr = 256, rmax = r_halo*4.d0, ns = 1024, a_ini = 1.d0/(1.d0+z_ini), delta_ini = delta_ini, r_halo = r_halo, bw_halo = r_halo/20.d0)
  call fig%open("rho.txt")
  call figphi%open("phi.txt")
  call fig%init(xlabel = "$r$", ylabel = "$\rho$")
  call figphi%init(xlabel = "$r$", ylabel = "$\phi$")
  call fig%plot(halo%r, halo%rho, color="blue", linetype="solid")
  call figphi%plot(halo%r, halo%phiave(1:halo%nr),color="blue", linetype="solid")
  dlna = log(1.d0+z_ini)/nsteps
  do while(halo%a  .lt. 0.5d0)
     call halo%evolve(dlna)
  enddo

  call fig%plot(halo%r, halo%rho, color="green", linetype="dashed", linewidth=1.5)
  call figphi%plot(halo%r, halo%phiave(1:halo%nr), color="green", linetype="dashed", linewidth=1.5)
  do while(halo%a .lt. 0.8d0)
     call halo%evolve(dlna)
  enddo
  call fig%plot(halo%r, halo%rho, color="red", linetype="dotted", linewidth=2.)
  call figphi%plot(halo%r, halo%phiave(1:halo%nr), color="red", linetype="dotted", linewidth=2.)
  call fig%close()
  call figphi%close()

end program fR1d
