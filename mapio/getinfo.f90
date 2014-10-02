program map
  use coop_healpix_mod
  use coop_wrapper_utils
  use healpix_types
  use alm_tools
  use pix_tools
  use head_fits
  use fitstools
  implicit none
#include "constants.h"
  
  integer(8) npixtot
  COOP_STRING::fin
  integer nside, nmaps, ordering, i, imap
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
        do imap = 1, nmaps
           call system("map2gif -inp "//trim(fin)//" -out "//trim(coop_file_name_of(fin, .false.))//"_"//COOP_STR_OF(imap)//".gif -bar T -sig "//COOP_STR_OF(imap))
        enddo
     else
        write(*,"(A)") trim(fin)//" does not exist"
     endif
     i = i+1
     fin = trim(coop_InputArgs(i))
  enddo
end program map
