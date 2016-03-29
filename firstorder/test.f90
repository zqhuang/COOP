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
  call coop_set_default_zeta_r(cosmology, cosmology%source(0))
end program test
