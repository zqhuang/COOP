program fR1d
  use coop_wrapper_firstorder
  use fR1d_mod
  implicit none
#include "constants.h"
  COOP_INT,parameter::nsave = 20, nr=3200
  COOP_REAL,parameter::r_halo = 0.001d0
  COOP_REAL,parameter::bw_halo = r_halo/40.d0, &
       rmax =  r_halo*20.d0, &
       dtau = rmax/nr/4.d0 
  COOP_REAL,parameter::  Omega_m = 1.d0,  &
       z_ini = 99.d0 , &
       delta_0 =  0.6d0*(1.5d0*coop_pi)**(2.d0/3.d0)*0.99
  type(coop_fr1d_obj)::halo
  COOP_INT::i, nsteps
  COOP_REAL::delta_ini, zsave(nsave), asave(nsave)
  COOP_STRING::prefix
  zsave = (/ 20.d0, 10.d0, 5.d0, 2.d0,  1.d0, 0.9d0, 0.8d0, 0.7d0, 0.6d0, 0.5d0, 0.45d0, 0.4d0, 0.35d0, 0.3d0, 0.25d0, 0.2d0, 0.15d0, 0.1d0, 0.05d0, 0.d0/)
  asave = 1.d0/(1.d0+zsave)
  write(*, "(A, E15.4)") "M = ",  0.5d0 * r_halo**3 * coop_SI_PlanckMass * coop_SI_hbyH0/ coop_SI_PlanckTime / coop_SI_Msun / 0.7d0
  delta_ini = delta_0 * coop_Growth_fitting(Omega_m, -1.d0, z_ini) / coop_Growth_fitting(Omega_m, -1.d0, 0.d0)
  delta_ini = delta_ini + 17.d0/21.d0*delta_ini**2 + ((17.d0/21.d0)**2*2.d0-341.d0/567.d0)*delta_ini**3
  call halo%init(Omega_m = Omega_m, nr = nr, rmax = rmax, a_ini = 1.d0/(1.d0+z_ini), delta_ini = delta_ini, r_halo = r_halo, bw_halo = bw_halo)
  if(trim(coop_InputArgs(1)).eq. "GR")then
     halo%do_GR = .true.
  else
     halo%do_GR = .false.
     if(trim(coop_InputArgs(1)).eq. "QS")then
        halo%QS_approx = .true.
     endif
  endif
  prefix = "col_"//trim(coop_InputArgs(1))//"_"
  do i=1, nsave
     do while(halo%a .lt. asave(i))
        call halo%evolve(dtau)
        if(mod(halo%nstep, 5000).eq.0)write(*, "(6G14.5)") halo%a, halo%lnrho(0),sign(log(abs(halo%phi(0:3))+exp(-100.d0)), halo%phi(0:3))
     enddo
     call halo%feedback(trim(prefix)//trim(coop_ndigits(i,2))//".txt")
  enddo
end program fR1d
