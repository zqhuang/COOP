program test
  use coop_wrapper_utils
  use coop_expint_mod  
  implicit none
#include "constants.h"

  type(coop_fits_file)::fp
  COOP_INT::nrows, ncols, index_minage, index_maxage, index_z, index_age, i
  COOP_REAL,dimension(:), allocatable::z, age, minage, maxage
  print*,coop_Threej000(2,3,1)**2*3.    
  print*,coop_Threej000(2,3,3)**2*7.
  print*,coop_Threej000(2,3,5)**2*11.  
  stop
  !  call fp%open('../../../lgal/portsmouth_stellarmass_passive_krou-DR12.fits')
  call fp%open('../../../lgal/portsmouth_stellarmass_passive_salp-DR12.fits')
  call fp%move_to_hdu(2)
  !  call fp%header%print()
  call fp%get_nrows_ncols(nrows, ncols)
  allocate(minage(nrows), maxage(nrows), z(nrows), age(nrows))
  index_minage = fp%get_col('MINAGE')
  index_maxage = fp%get_col('MAXAGE')
  index_z = fp%get_col('Z')
  index_age = fp%get_col('AGE')
  call fp%load_double_column(index_minage, minage)
  call fp%load_double_column(index_maxage, maxage)
  call fp%load_double_column(index_z, z)
  call fp%load_double_column(index_age, age)
  call fp%close()
  do i=1, nrows
     write(*,'(4F11.4)')  z(i), age(i), maxage(i)-age(i),  age(i) - minage(i)
  enddo
!  print*, sqrt(sum((maxage-age)**2)/nrows),  sqrt(sum((age-minage)**2)/nrows)
end program test
