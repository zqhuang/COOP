program fR1d
  use coop_wrapper_firstorder
  use fR1d_mod
  implicit none
#include "constants.h"
  COOP_INT,parameter::nsave = 5
  type(coop_fr1d_obj)::halo
  type(coop_file)::frho, fphi
  COOP_INT::i, nsteps
  COOP_REAL::delta_ini, z_ini, delta_0, Omega_m, r_halo, zsave(nsave), asave(nsave)
  COOP_REAL,dimension(:,:),allocatable::lnrho, lnphi


  Omega_m = 1.d0
  z_ini = 99.d0
  delta_0 =  0.6d0*(1.5d0*coop_pi)**(2.d0/3.d0)*0.8
  print*, "delta_0 = ", delta_0
  r_halo = 0.001d0
  zsave = (/ 98.d0, 50.d0, 10.d0, 2.d0, 1.d0 /)
  asave = 1.d0/(1.d0+zsave)
  write(*, "(A, E15.4)") "M = ",  0.5d0 * r_halo**3 * coop_SI_PlanckMass * coop_SI_hbyH0/ coop_SI_PlanckTime / coop_SI_Msun / 0.7d0
  delta_ini = delta_0 * coop_Growth_fitting(Omega_m, -1.d0, z_ini) / coop_Growth_fitting(Omega_m, -1.d0, 0.d0)
  delta_ini = delta_ini + 17.d0/21.d0*delta_ini**2
  call halo%init(Omega_m = Omega_m, nr = 1024, rmax = r_halo*20.d0, a_ini = 1.d0/(1.d0+z_ini), delta_ini = delta_ini, r_halo = r_halo, bw_halo = r_halo/10.d0, dtau = 2.d-5)
  halo%ignore_cameleon_force = .false.
  halo%QS_approx = .true.

  allocate(lnrho(0:halo%nr, nsave), lnphi(0:halo%nr, nsave))
  call frho%open("rho.txt")
  call fphi%open("phi.txt")
  do i=1, nsave
     do while(halo%a .lt. asave(i))
        call halo%evolve()
     enddo
     lnrho(:, i) = halo%lnrho(:, halo%ind_u(1))
     lnphi(:, i) = halo%lnphi(:, halo%ind_phi(1))
  enddo
  do i=0, halo%nr
     write(frho%unit, "("//COOP_STR_OF(nsave+1)//"E16.7)") halo%r(i), lnrho(i, :)
     write(fphi%unit, "("//COOP_STR_OF(nsave+1)//"E16.7)") halo%r(i), lnphi(i, :)
  enddo
  call frho%close()
  call fphi%close()
end program fR1d
