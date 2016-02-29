program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  COOP_INT::numhdus, ihdu
  type(coop_dictionary)::header
  COOP_STRING::filename
  COOP_INT,parameter::n9=2479
  COOP_SHORT_INT::ells(2508)
  COOP_SINGLE::Dls(2508), Els(2508), Dls_theory(2508), Eb, weight
  COOP_INT::l, l_start, l_end
  type(coop_file)::fp
  call fp%open("planck14best_lensedCls.dat")
  do l=2, 2500
     read(fp%unit, *) ells(l), Dls_theory(l)
  enddo
  filename = "planck14/COM_PowerSpect_CMB_R2.02.fits"
  call coop_convert_to_C_string(filename)
  call coop_fits_read_col_short(filename, 9, 1, n9, ells(30:))
  call coop_fits_read_col_float(filename, 9, 2, n9, Dls(30:))
  call coop_fits_read_col_float(filename, 9, 3, n9, Els(30:))
  write(*,*) "enter l_start, l_end:"
30  read(*,*) l_start, l_end
  if(l_end .le. l_start) stop
  Eb = 0.d0
  weight = 0.d0
  do l = l_start, l_end
     Eb = Eb + (Dls(l) - Dls_theory(l))/Els(l)**2
     weight = weight + 1.d0/Els(l)**2
  enddo
  Eb = Eb/weight
  write(*,*) (Eb*sqrt(weight))**2
  goto 30
end program test
