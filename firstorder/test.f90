program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::cosmology
  type(coop_real_table)::paramtable
  logical success
  COOP_INT,parameter::nk = 300
  COOP_REAL::k(nk), kMpc(nk), Dz(nk), z
  type(coop_asy)::fig
  COOP_INT::i
  call paramtable%insert("ombh2", 0.01d0)
  call paramtable%insert("omch2", 0.45d0)
  call paramtable%insert("h", 0.7d0)
  call paramtable%insert("tau", 0.08d0)
  call paramtable%insert("As", 2.1d-9)
  call paramtable%insert("ns", 0.967d0)
  call cosmology%set_up(paramtable, success)
  if(.not. success) stop "initialization failed"
  call cosmology%compute_source(0, success)
  call fig%open("growth.txt")
  call fig%init(xlabel = "$k\,[\mathrm{Mpc}^{-1}]$", ylabel = "$(1+z)D(z)|_{z=0}$", xlog = .true.)
  call coop_set_uniform(nk, kMpc, 0.0001d0, 0.2d0, logscale = .true.)
  k = kMpc/cosmology%H0Mpc()
  z = 25.d0
  do i = 1, nk
     Dz(i) = (cosmology%growth_of_z(k = k(i), z = z)*(1.d0+z)) 
  enddo
  call fig%plot(k, Dz)
  call fig%close()
end program test
