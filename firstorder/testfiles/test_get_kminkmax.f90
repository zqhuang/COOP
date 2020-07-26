program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::cosmology
  type(coop_real_table)::paramtable
  logical success
  COOP_REAL::fsky
  COOP_INT::nz, i
  COOP_REAL:: z(100), dz(100), k_min(100), k_max(100)
  call paramtable%insert("ombh2", 0.022d0)
  call paramtable%insert("omch2", 0.12d0)
  call paramtable%insert("h", 0.67d0)
  call paramtable%insert("tau", 0.07d0)
  call paramtable%insert("As", 2.15d-9)
  call paramtable%insert("ns", 0.969d0)
  call cosmology%set_up(paramtable, success)
  if(.not. success) stop "initialization failed"
  call cosmology%compute_source(0, success)
  write(*,*) "enter fsky:"
  read(*,*) fsky
  write(*,*) "number of redshifts:"
  read(*,*) nz
  write(*,*) "enter z list"
  read(*,*) z(1:nz)
  write(*,*) "enter dz list"
  read(*, *) dz(1:nz)
  do i = 1, nz
     k_max(i) = kmax(z(i))
     k_min(i) = kmin(z(i), dz(i))
  enddo
  write(*,'(A, '//COOP_STR_OF(nz)//'G15.6)') 'kmin = ', k_min(1:nz)
  write(*,'(A, '//COOP_STR_OF(nz)//'F10.6)') 'kmax = ', k_max(1:nz)
  do i = 1, nz
     write(*,"(A, G15.6)") "window"//COOP_STR_OF(i)//" = ", sqrt(coop_ln2)/coop_pi*k_min(i)
  enddo

contains

  function Volume(z, dz) result(V)
    COOP_REAL::z, dz, V
    V = coop_4pi/3.d0*fsky* (cosmology%comoving_angular_diameter_distance(1.d0/(1.d0+z+dz/2.d0))**3 - cosmology%comoving_angular_diameter_distance(1.d0/(1.d0+z-dz/2.d0))**3) * (3000.)**3
  end function Volume

  function numz(z, dz) 
    COOP_REAL::z, dz, numz
    numz = coop_integrate(dndz, z-dz/2.d0, z+dz/2.d0)*(coop_4pi*fsky/coop_SI_arcmin**2) / Volume(z, dz)
  end function numz

  function dndz(z) 
    COOP_REAL::z, dndz
    dndz = 640 *z**2*exp(-z/0.35)
  end function dndz

  function kmin(z, dz)
    COOP_REAL::z, kmin, dz
    kmin = coop_2pi/(Volume(z, dz)*(3.d0/coop_4pi))**(1.d0/3.d0)
  end function kmin

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
       if(sigmamid .gt. 0.501d0)then 
          rlow =rmid
       elseif(sigmamid .lt. 0.499d0)then
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
