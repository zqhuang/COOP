program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::cosmology
  COOP_REAL::a, b
  type(coop_real_table)::paramtable
  logical success
  COOP_INT, parameter::nw = 100
  COOP_REAL,parameter::sigma_W = 5.d0
  COOP_INT i
  COOP_REAL::k, kw(nw), wsq(nw), z, pk1, pk2, pk3
  call paramtable%insert("ombh2", 0.022d0)
  call paramtable%insert("omch2", 0.11d0)
  call paramtable%insert("theta", 1.04d0)
  call paramtable%insert("tau", 0.07d0)
  call paramtable%insert("As", 2.d-9)
  call paramtable%insert("ns", 0.96d0)
  call coop_prtsystime(.true.)
  write(*,*) "set up"
  call cosmology%set_up(paramtable, success)
  if(.not. success) stop "initialization failed"
  print*, cosmology%cosmomc_theta()
  call coop_prtsystime()
  write(*,*) "solve pert"
  call coop_prtsystime(.true.)
  call cosmology%compute_source(0, success)
  if(.not. success) stop "linear pert failed"
  call coop_prtsystime()
!!$
!!$  write(*,*) "compute cls"
!!$  call coop_prtsystime(.true.)
!!$  call cosmology%set_cls(0, 2, 3000)
!!$  call coop_prtsystime()
!!$  a=cosmology%source(0)%Cls(coop_index_ClTT, 200)
!!$  write(*,*)"update cls"
!!$  call coop_prtsystime(.true.)
!!$  call paramtable%update("As", 2.2d-9)
!!$  call cosmology%set_primordial_power(paramtable)
!!$  call cosmology%update_cls(0)
!!$  call coop_prtsystime()
  do i=1, nw
     kw(i) = (i-0.99d0)* (5.d0*sigma_W/nw)
     wsq(i) = (((kw(i)/sigma_W)*cos(kw(i)/sigma_W) - sin(kw(i)/sigma_W))/(kw(i)/sigma_W)**3)**2
  enddo
  z = 0.5d0
  do i = 1, 50

     k = sigma_W*exp((i-1)*0.2d0)
     pk1 = cosmology%Gaussian_smeared_matter_power(z= z, k = k, sigma_W = sigma_W)
     pk2 = cosmology%smeared_matter_power(z = z, k = k, nw = nw, kw = kw, wsq = wsq)
     pk3 = cosmology%matter_power(z = z, k = k)
     write(*,*) k, pk1, pk2, pk3
  enddo
end program test
