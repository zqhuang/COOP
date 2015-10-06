program test
  use coop_wrapper_firstorder
  use coop_halofit_mod
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::cosmology
  type(coop_file)fp
  COOP_INT, parameter::nk=300, nz = 200, num_l = 150, nbins = 1
  COOP_REAL,parameter::alpha = 2.d0, beta = 1.5d0, z0= 0.667d0
  COOP_INT ik
  COOP_REAL k(nk), Pl(nk), Pnl(nk), Pl_weyl(nk), Pnl_weyl(nk), zmean
  COOP_REAL z_source(nz), n_source(nz, nbins), l(num_l), Cl(num_l, nbins*(nbins+1)/2)

  !!set cosmology
  call cosmology%Set_Planck_bestfit()
  call cosmology%compute_source(0)

!!$  call coop_set_uniform(nk, k, 0.001d0, 1.d4, logscale = .true.)
!!$  k = k/cosmology%H0Mpc()
!!$
!!$  call coop_halofit_get_power(cosmology, 0.d0, nk, k, Pl, Pnl)
!!$  call coop_halofit_get_Weyl_power(cosmology, 0.d0, nk, k, Pl_Weyl, Pnl_weyl)
!!$  k = k * cosmology%H0Mpc()  !!
!!$  Pl = Pl/k**3*(2.*coop_pi2)
!!$  Pnl = Pnl/k**3*(2.*coop_pi2)
!!$  
!!$  Pl_Weyl = Pl_Weyl/k**3*(2.*coop_pi2)
!!$  Pnl_Weyl = Pnl_Weyl/k**3*(2.*coop_pi2)
!!$
!!$  call fp%open("pnl.txt", "w")
!!$  do ik=1,nk
!!$     write(fp%unit, "(20E16.7)") k(ik), Pl(ik), Pnl(ik), Pl_weyl(ik), Pnl_weyl(ik)
!!$  enddo
!!$  call fp%close()
!!$stop
  zmean = z0*gamma((2.d0+alpha)/beta)/gamma((1.d0+alpha)/beta)
  print*, "<z> = ", zmean
  call coop_set_uniform(nz, z_source, 0.01d0, zmean*3.d0)
  call coop_set_uniform(num_l, l, 1.d0, 1.d4, logscale = .true.)
  n_source(:, 1) =(z_source/z0)**alpha*exp(-(z_source/z0)**beta)/(z0/beta*Gamma((1.d0+alpha)/beta))
  print*, sum(n_source)*(z_source(2)-z_source(1))
  call coop_halofit_get_WeakLensing_Power(cosmology, nbins, nz, z_source, n_source, num_l, l, Cl)

  call fp%open("WLCl.txt", "w")
  do ik=1,num_l
     write(fp%unit, "(20E16.7)") l(ik), l(ik)**2*Cl(ik,:)/coop_2pi
  enddo
  call fp%close()
  
  
end program test
