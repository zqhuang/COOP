program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::cosmology
  type(coop_real_table)::paramtable
  logical success
  COOP_REAL::fsky = 0.558
  COOP_REAL:: z, dz, k
  call paramtable%insert("ombh2", 0.022d0)
  call paramtable%insert("omch2", 0.12d0)
  call paramtable%insert("theta", 1.04d0)
  call paramtable%insert("tau", 0.08d0)
  call paramtable%insert("As", 2.1d-9)
  call paramtable%insert("ns", 0.967d0)
  call cosmology%set_up(paramtable, success)
  if(.not. success) stop "initialization failed"
  call cosmology%compute_source(0, success)
  do 
     read(*, *) z
     k = kmax(z)
     write(*,*) k
  enddo

contains

  function Volume(z, dz) result(V)
    COOP_REAL::z, dz, V
    V = coop_4pi/3.d0*fsky* (cosmology%comoving_angular_diameter_distance(1.d0/(1.d0+z+dz/2.d0))**3 - cosmology%comoving_angular_diameter_distance(1.d0/(1.d0+z-dz/2.d0))**3) * (3000.)**3
  end function Volume

  function nz(z, dz) 
    COOP_REAL::z, dz, nz
    nz = coop_integrate(dndz, z-dz/2.d0, z+dz/2.d0)*(coop_4pi*fsky/coop_SI_arcmin**2) / Volume(z, dz)
  end function nz

  function dndz(z) 
    COOP_REAL::z, dndz
    dndz = 640 *z**2*exp(-z/0.35)
  end function dndz

  function kmax(z)
    COOP_REAL::z, kmax, rlow, rhigh, sigmalow, sigmahigh, sigmamid, rmid
    rlow = 0.5d0
    rhigh = 20.d0
    sigmalow = cosmology%sigma_tophat_R(z=z, r = rlow/cosmology%h()*cosmology%H0Mpc())
    sigmahigh = cosmology%sigma_tophat_R(z=z, r = rhigh/cosmology%h()*cosmology%H0Mpc())
    if(sigmalow .lt. 0.5d0 .or. sigmahigh .gt. 0.5d0)stop "kmax failed"
    do while(rhigh/rlow .gt. 1.002)
       rmid = sqrt(rlow*rhigh)
       sigmamid = cosmology%sigma_tophat_R(z= z, r  = rmid/cosmology%h()*cosmology%H0Mpc())
       if(sigmamid .gt. 0.5001d0)then 
          rlow =rmid
       elseif(sigmamid .lt. 0.4999d0)then
          rhigh = rmid
       else
          kmax = 1.d0/rmid
          return
       endif
    enddo
    rmid = sqrt(rlow*rhigh)
    kmax = coop_pio2/rmid
  end function kmax

end program test
