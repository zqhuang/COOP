program shells
  use coop_hnn_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask
  integer i
  COOP_INT::fsky 
  COOP_REAL::cut_lower, cut_upper, cut_mid, dcut, fsky_mid
  write(*,*) "Enter fsky = "
  read(*,*) fsky
  if(fsky .eq. 80)then  !!initialization
     call map%read("plmaps/CO_5deg.fits", nmaps_wanted = 1)
     call mask%init(nside= 2048, nmaps=1, genre="MASK")
     mask%map = 1.
     where(map%map .gt. 4.e5)
        mask%map = 0.
     end where
     call mask%write("plmaps/mask_LR80_n2048.fits")
  else
     call map%read("plmaps/I857_5deg.fits", nmaps_wanted=1)
     call mask%read("plmaps/mask_LR80_n2048.fits", nmaps_wanted=1)
     print*,sum(dble(mask%map))/dble(mask%npix)
     cut_lower = minval(map%map)
     cut_upper = maxval(map%map)
     dcut = (cut_upper - cut_lower) * 1.e-5
     where(mask%map .eq. 0.)
        map%map = cut_upper+1.
     end where
     do while(cut_upper - cut_lower .gt. dcut)
        cut_mid = (cut_lower + cut_upper)/2.
        fsky_mid = count( map%map .lt. cut_mid )/dble(map%npix)
        if(fsky_mid .gt. fsky/100.)then
           cut_upper = cut_mid
        else
           cut_lower = cut_mid
        endif
        print*, cut_upper, cut_lower, cut_mid, fsky_mid
     enddo
     cut_mid = (cut_lower + cut_upper)/2.
     where(map%map(:,1) .gt. cut_mid)
        mask%map(:,1) = 0.
     end where
     print*,sum(dble(mask%map))/dble(mask%npix)
     call mask%write("plmaps/mask_LR"//COOP_STR_OF(fsky)//"_n2048.fits")
  endif
end program shells
