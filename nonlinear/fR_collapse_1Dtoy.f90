program fR1d
  use coop_wrapper_firstorder
  use fR1d_mod
  implicit none
#include "constants.h"
  COOP_INT,parameter::nsave = 7
  type(coop_fr1d_obj)::halo
  COOP_INT::i, nsteps
  COOP_REAL::delta_ini, z_ini, delta_0, Omega_m, r_halo, zsave(nsave), asave(nsave)
  COOP_REAL,dimension(:,:),allocatable::lnrho, lnphi

  Omega_m = 1.d0
  z_ini = 99.d0
  delta_0 =  0.6d0*(1.5d0*coop_pi)**(2.d0/3.d0)
  print*, "delta_0 = ", delta_0
  r_halo = 0.001d0
  zsave = (/ 1.d0, 0.5d0, 0.4d0, 0.3d0, 0.2d0, 0.1d0, 0.d0/)
  asave = 1.d0/(1.d0+zsave)
  write(*, "(A, E15.4)") "M = ",  0.5d0 * r_halo**3 * coop_SI_PlanckMass * coop_SI_hbyH0/ coop_SI_PlanckTime / coop_SI_Msun / 0.7d0
  delta_ini = delta_0 * coop_Growth_fitting(Omega_m, -1.d0, z_ini) / coop_Growth_fitting(Omega_m, -1.d0, 0.d0)
  delta_ini = delta_ini + 17.d0/21.d0*delta_ini**2
  call halo%init(Omega_m = Omega_m, nr = 1024, rmax = r_halo*20.d0, a_ini = 1.d0/(1.d0+z_ini), delta_ini = delta_ini, r_halo = r_halo, bw_halo = r_halo/20.d0, dtau = r_halo/200.d0)
  halo%QS_approx = .true.
  halo%do_GR = .false.
  allocate(lnrho(0:halo%nr, nsave), lnphi(0:halo%nr, nsave))
  do i=1, nsave
     do while(halo%a(halo%time(1)) .lt. asave(i))
        call halo%evolve()
        if(mod(halo%nstep, 5000).eq.0)print*, halo%a(halo%time(1)), halo%lnrho(1, halo%time(1))
     enddo
     call halo%feedback("save_QS"//COOP_STR_OF(i)//".txt")
     lnrho(:, i) = halo%lnrho(:, halo%time(1))
     if(.not.halo%do_GR)lnphi(:, i) = halo%lnphi(:, halo%time(1))
  enddo
end program fR1d
