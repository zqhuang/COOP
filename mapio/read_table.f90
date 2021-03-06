program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  use coop_fitsio_mod
  implicit none
#include "constants.h"
  COOP_INT::numhdus, ihdu
  type(coop_dictionary)::header
  COOP_STRING::filename
  COOP_INT,parameter::n9=2479
  COOP_INT::ells(2508)
  COOP_REAL::Dls(2508), Els(2508), Dls_theory(2508), Eb, weight
  COOP_INT::l, l_start, l_end, i
  type(coop_file)::fp
  type(coop_fits_file)::ff
  call fp%open("planck14best_lensedCls.dat")
  do l=2, 2500
     read(fp%unit, *) ells(l), Dls_theory(l)
  enddo
  call fp%close()
  call ff%open("planckCls.fits", ihdu=1) !"planck14/COM_PowerSpect_CMB_R2.02.fits", ihdu = 9)
  call ff%load_int_column(col = 1, data = ells(30:))
  call ff%load_double_column(col = 2, data = Dls(30:))
  call ff%load_double_column(col = 3, data = Els(30:))
  call ff%close()

  write(*,*) "enter l_start, l_end:"
30  read(*,*) l_start, l_end
  if(l_end .le. l_start .or. l_start .lt. 30) stop
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
