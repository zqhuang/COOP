program test
  implicit none
#include "constants.h"
  integer nx, ny, bytes
  real(8),dimension(:,:),allocatable::data
  character(LEN=*),parameter::mapdir = "../data/cmb/maps/act/"
  character(LEN=*),parameter::fitsfile = mapdir//"Q50.fits"
  call coop_fits_print_header(fitsfile)
  call coop_fits_get_dimension(fitsfile, nx, ny, bytes)
  print*, nx, ny, bytes
  if(bytes.eq.8)then
     allocate(data(ny, nx))
     call coop_fits_get_double_data(fitsfile, data, nx, ny)
     print*, data(100, 100)
  endif
end program test
