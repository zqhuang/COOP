program read_table
  implicit none
  character(LEN=1024)::file_root, inifile
  real*8, allocatable::zcol(:,:,:)
  real*8 fpk_min, fpk_max, e_min, e_max, pbye_min, pbye_max, zmax
  integer n_fpk, n_e, n_pbye, nz, iz
  real*8, dimension(:),allocatable::z, t, H, chi, D, H_D

  file_root  = "myoutputs/test"
  open(10, FILE = trim(adjustl(file_root))//"_zcol.dat", access = 'stream')
  read(10) n_fpk, fpk_min, fpk_max
  read(10) n_e, e_min, e_max
  read(10) n_pbye, pbye_min, pbye_max
  allocate(zcol(n_fpk, n_e, n_pbye))
  read(10) zcol
  close(10)
  open(10, FILE = trim(adjustl(file_root))//"_background.dat", access = 'stream')
  read(10) nz, zmax
  allocate(z(nz), t(nz), H(nz), chi(nz), D(nz), H_D(nz))
  do iz=1, nz
     read(10) z(iz), t(iz), H(iz), chi(iz), D(iz), H_D(iz)
  enddo
  close(10)

!!units: 
!!t : 1/H_0
!!H : H_0
!!chi:  c/H_0 (comoving distance)
!!D :   normalized to 1 at z = 0
!!H_D:  H_0 (defined as dln D/ dt)

end program read_table
