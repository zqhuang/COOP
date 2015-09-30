program map
  use coop_healpix_mod
  use coop_wrapper_utils
#ifdef HAS_HEALPIX  
  use healpix_types
  use alm_tools
  use pix_tools
  use head_fits
  use fitstools
#endif  
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX  
  integer(8) npixtot
  COOP_STRING::fin, gif
  COOP_REAL mean
  integer nside, nmaps, ordering, i, imap
  type(coop_healpix_maps)::hgm
  i = 1
  fin = trim(coop_InputArgs(i))
  do while(trim(fin).ne."")
     if(coop_file_exists(fin))then
        npixtot = getsize_fits(trim(fin), nmaps = nmaps, nside = nside, ordering = ordering)
        write(*,"(A)")"**********************************"
        write(*, "(A)")  "file: "//trim(fin)
        write(*, "(A,I5)") "nmaps = "//COOP_STR_OF(nmaps)
        write(*, "(A,I5)") "nside = "//COOP_STR_OF(nside)
        if(ordering .eq. COOP_RING)then
           write(*, "(A)") "ordering: ring"
        elseif(ordering .eq. COOP_NESTED)then
           write(*, "(A)") "ordering: nested"
        else
           write(*, "(A)") "ordering: unknown"
        endif
        call hgm%read(trim(fin))
        do imap = 1, nmaps
           mean = sum(hgm%map(:, imap))/hgm%npix
           write(*, "(A)") "--- map #"//COOP_STR_OF(imap)//" --- "
           write(*, "(A, I5)") "spin = ", hgm%spin(imap)
           write(*, "(A)") "    #pix>0 fsky= "//COOP_STR_OF(count(hgm%map(:, imap).gt. 0.)/dble(hgm%npix))
           write(*, "(A)") "    #pix<0 fsky= "//COOP_STR_OF(count(hgm%map(:, imap).lt. 0.)/dble(hgm%npix))
           write(*, "(A)") "    min = "//COOP_STR_OF(minval(hgm%map(:,imap)))
           write(*, "(A)") "    max = "//COOP_STR_OF(maxval(hgm%map(:,imap)))           
           write(*, "(A)") "    mean = "//COOP_STR_OF(mean)
           write(*, "(A)") "    rms = "//COOP_STR_OF(sqrt(sum((hgm%map(:, imap)-mean)**2)/hgm%npix))
!!$           gif = "tempgifs/"//trim(coop_file_name_of(fin, .false.))//"_"//COOP_STR_OF(imap)//".gif"
!!$           if(.not. coop_file_exists(gif)) &
!!$                call system("map2gif -inp "//trim(fin)//" -out "//trim(gif)//" -bar T -sig "//COOP_STR_OF(imap))
        enddo
        call hgm%free()
     else
        write(*,"(A)") trim(fin)//" does not exist"
     endif
     i = i+1
     fin = trim(coop_InputArgs(i))
  enddo
#else
  stop "you need to install healpix"
#endif
end program map
