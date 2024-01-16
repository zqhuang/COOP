program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map
  type(coop_file)::fp
  COOP_INT::i, j, l
  COOP_INT,parameter::lmin = 19, lmax = 31
  COOP_REAL::Cls_model(2:lmax), Cls(2:lmax), junk, Cl_binned, weights(2:lmax)
  COOP_STRING::line
  do l = lmin, lmax
     weights(l) = (2.d0*l+1.d0)
  enddo
  call map%init(nmaps = 1, nside = 2, genre = "T", nested = .true.)
  do i=0, map%npix-1
     call fp%open("jycls/inpainted_cls_ipix"//COOP_STR_OF(i)//".dat", "r")
     read(fp%unit, "(A)") line
     do j= 2, 31
        read(fp%unit, *) l, Cls_model(l), Cls(l), junk
     enddo
     Cl_binned = sum(Cls(lmin:lmax)*weights(lmin:lmax))/sum(Cls_model(lmin:lmax)*weights(lmin:lmax)) - 1.d0
     call fp%close()
     map%map(i, 1) = Cl_binned
  enddo
  call map%write("deltaCl.fits")
end program test
