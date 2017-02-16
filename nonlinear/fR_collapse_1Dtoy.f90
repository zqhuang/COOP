program fR1d
  use coop_wrapper_firstorder
  use fR1d_mod
  implicit none
#include "constants.h"
  COOP_INT,parameter::nsave = 5
  type(coop_fr1d_obj)::halo, haloQS
  type(coop_file)::frho, fphi
  COOP_INT::i, nsteps
  COOP_REAL::delta_ini, z_ini, delta_0, Omega_m, dlna, r_halo, zsave(nsave), asave(nsave)
  COOP_REAL,dimension(:,:),allocatable::rho, phi

  Omega_m = 0.3d0
  z_ini = 99.d0
  delta_0 = 2.d0
  r_halo = 0.001d0
  zsave = (/ 10.d0, 2.d0, 1.d0, 0.5d0, 0.d0 /) 
  asave = 1.d0/(1.d0+zsave)
  write(*, "(A, E15.4)") "M = ",  0.5d0 * r_halo**3 * coop_SI_PlanckMass * coop_SI_hbyH0/ coop_SI_PlanckTime / coop_SI_Msun / 0.7d0
  delta_ini = delta_0 * coop_Growth_fitting(Omega_m, -1.d0, z_ini) / coop_Growth_fitting(Omega_m, -1.d0, 0.d0)
  call halo%init(Omega_m = Omega_m, nr = 512, rmax = r_halo*5.d0, ns = 2048, a_ini = 1.d0/(1.d0+z_ini), delta_ini = delta_ini, r_halo = r_halo, bw_halo = r_halo/20.d0)
  dlna = -log(halo%a)/20000
  haloQS = halo
  haloQS%QS_approx = .true.
  halo%QS_approx = .false.
  allocate(rho(halo%nr, 2*nsave), phi(halo%nr, 2*nsave))
  call frho%open("rho.txt")
  call fphi%open("phi.txt")
  do i=1, nsave
     do while(halo%a .lt. asave(i))
        call halo%evolve(dlna)
     enddo
     do while(haloQS%a .lt. asave(i))
        call haloQS%evolve(dlna)
     enddo
     rho(:, 2*i-1) = halo%rho(1:halo%nr)
     phi(:, 2*i-1) = halo%phi(1:halo%nr, halo%time(1))
     rho(:, 2*i) = haloQS%rho(1:halo%nr)
     phi(:, 2*i) = haloQS%phi(1:halo%nr, haloQS%time(1))
  enddo
  do i=1, halo%nr
     write(frho%unit, "("//COOP_STR_OF(2*nsave+1)//"E16.7)") halo%r(i), rho(i, :)
     write(fphi%unit, "("//COOP_STR_OF(2*nsave+1)//"E16.7)") halo%r(i), phi(i, :)
  enddo
  call frho%close()
  call fphi%close()
end program fR1d
