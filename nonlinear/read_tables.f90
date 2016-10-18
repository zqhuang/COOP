program test
  implicit none
  character(LEN=256)::fname
  real*8 fpk_min, fpk_max, e_min, e_max, pbye_min, pbye_max, dlnfpk, de, dpbye, input_fpk, input_e, input_p
  integer n_fpk, n_e, n_pbye
  real*8,allocatable::zcol(:,:,:)

  real*8::H0, Omega_b, Omega_c, A_s, n_s, sigma_8, w, wa
  real*8 zmax, dz, input_z
  integer nz, iz
  real*8, dimension(:),allocatable::z, t, H, chi, D, H_D

  !!load z_collapse table
  write(*,*) "Enter the file name of the z_collapse table:"
  read(*,*) fname
  open(10, file=fname, access ='stream')
  read(10) n_fpk, fpk_min, fpk_max
  read(10) n_e, e_min, e_max
  read(10) n_pbye, pbye_min, pbye_max
  allocate(zcol(n_fpk, n_e, n_pbye))
  read(10) zcol
  close(10)

  !!load background table
  write(*,*) "Enter the file name of the background evolution table:"
  read fname
  open(20, file=fname, access ='stream')
  read(20)H0, Omega_b, Omega_c, A_s, n_s, sigma_8, w, wa
  read(20) nz, zmax
  allocate(z(nz), t(nz), H(nz), chi(nz), D(nz), H_D(nz))
  do iz = 1, nz
     read(20) z(iz), t(iz), H(iz), chi(iz), D(iz), H_D(iz)
  enddo
  close(20)

  dlnfpk = log(fpk_max/fpk_min)/(n_fpk-1.d0)
  de = (e_max - e_min)/(n_e-1.d0)
  dpbye = (pbye_max - pbye_min)/(n_pbye-1.d0)
  dz = zmax/(nz-1.d0)

  write(*,*) "Now test the z_collapse table."
  write(*,*) "+++++ Enter f_pk, e, p (enter an f_pk < 0 to stop) ++++++++++"
  read(*,*) input_fpk, input_e, input_p
  do while(input_fpk .gt. 0.d0)
     write(*,*) "z_collapse = ", zcol_interpolate(input_fpk, input_e, input_p)
     write(*,*) "+++++ Enter f_pk, e, p (enter an f_pk < 0 to stop) ++++++++++"
     read(*,*) input_fpk, input_e, input_p
  enddo

  write(*,*) "Now test the background table."
  write(*,*) "z_max =", zmax
  write(*,*) "+++++++++ Enter the redshift (enter z<0 or z>=zmax to stop) +++++++++++"
  read(*,*) input_z
  do while(input_z .ge. 0.d0 .and. input_z .lt. zmax)
     write(*,*) "H / H_0 = ", H_interpolate(input_z)
     write(*,*) "chi / (c/H_0) = ", chi_interpolate(input_z)
     write(*,*) "t / H_0^{-1}  = ", t_interpolate(input_z)
     write(*,*) "D = ", D_interpolate(input_z)
     write(*,*) "H_D / H_0 = ", H_D_interpolate(input_z)
     write(*,*) "+++++++++ Enter the redshift (enter z<0 or z>=zmax to stop) +++++++++++"
     read(*,*) input_z
  enddo

contains

  function zcol_interpolate(fpk, e, p) result(z)
    real*8::fpk, e, p, z
    real*8::pbye,  rfpk, re, rp
    integer::ifpk, ie, ip
    if(e.ne.0.d0)then
       pbye = p/e
    else
       if(p.ne.0.)then
          z = -1.
          return
       endif
       pbye = 0.
    endif
    rfpk = log(fpk/fpk_min)/dlnfpk + 1.d0
    ifpk = floor(rfpk)
    if(ifpk .lt. 1 .or. ifpk .ge. n_fpk)then
       z = -1
       return
    endif
    rfpk = rfpk - ifpk
    re = (e-e_min)/de + 1.d0
    ie = floor(re)
    if(ie .lt. 1 .or. ie .ge. n_e)then
       z = -1
       return
    endif
    re = re - ie
    rp = (pbye - pbye_min)/dpbye + 1.d0
    ip = floor(rp) 
    if(ip .lt. 1 .or. ip .ge. n_pbye)then
       z = -1
       return
    endif
    rp = rp - ip

    z = (1.d0-rfpk)*((1.d0-re)*((1.d0-rp)*zcol(ifpk, ie, ip) + rp*zcol(ifpk, ie, ip+1) ) &
         +re*((1.d0-rp)*zcol(ifpk, ie+1, ip) + rp*zcol(ifpk, ie+1, ip+1) )) &
         +rfpk*((1.d0-re)*((1.d0-rp)*zcol(ifpk+1, ie, ip) + rp*zcol(ifpk+1, ie, ip+1) ) &
         +re*((1.d0-rp)*zcol(ifpk+1, ie+1, ip) + rp*zcol(ifpk+1, ie+1, ip+1) ))
  end function zcol_interpolate

  function H_interpolate(z) result(Hofz)
    real*8 z, Hofz, rz
    integer iz
    rz = z/dz+1.d0
    iz = floor(rz)
    if(iz .gt. nz)then
       Hofz = -1
       return
    endif
    rz = rz - iz
    Hofz = H(iz)*(1.d0-rz)+H(iz+1)*rz
  end function H_interpolate

  function chi_interpolate(z) result(chiofz)
    real*8 z, chiofz, rz
    integer iz
    rz = z/dz+1.d0
    iz = floor(rz)
    if(iz .gt. nz)then
       chiofz = -1
       return
    endif
    rz = rz - iz
    chiofz = chi(iz)*(1.d0-rz)+chi(iz+1)*rz
  end function chi_interpolate

  function H_D_interpolate(z) result(H_Dofz)
    real*8 z, H_Dofz, rz
    integer iz
    rz = z/dz+1.d0
    iz = floor(rz)
    if(iz .gt. nz)then
       H_Dofz = -1
       return
    endif
    rz = rz - iz
    H_Dofz = H_D(iz)*(1.d0-rz)+H_D(iz+1)*rz
  end function H_D_interpolate

  function t_interpolate(z) result(tofz)
    real*8 z, tofz, rz
    integer iz
    rz = z/dz+1.d0
    iz = floor(rz)
    if(iz .gt. nz)then
       tofz = -1
       return
    endif
    rz = rz - iz
    tofz = t(iz)*(1.d0-rz)+t(iz+1)*rz
  end function t_interpolate

  function D_interpolate(z) result(Dofz)
    real*8 z, Dofz, rz
    integer iz
    rz = z/dz+1.d0
    iz = floor(rz)
    if(iz .gt. nz)then
       Dofz = -1
       return
    endif
    rz = rz - iz
    Dofz = D(iz)*(1.d0-rz)+D(iz+1)*rz
  end function D_interpolate



end program test
