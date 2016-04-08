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
  COOP_STRING::fin, gif, fmask
  COOP_REAL mean, rms, mapmax, mapmin
  integer nside, nmaps, ordering, i, imap
  character::yn
  type(coop_healpix_maps)::hgm, mask
  logical::has_mask
  COOP_REAL::weight
  i = 1
  fin = trim(coop_InputArgs(i))
  i = 2
  fmask = trim(coop_InputArgs(i))
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
     if(trim(fmask).ne."")then
        if(coop_file_exists(fmask))then
           has_mask = .true.
           call mask%read(fmask)
        else
           write(*,*) "Error: mask file "//trim(fmask)//" does not exist."
           stop
        endif
     else
        has_mask = .false.
     endif
     do imap = 1, nmaps
        if(has_mask)then
           hgm%map(:,imap) = hgm%map(:, imap)*mask%map(:, 1)
           weight = sum(dble(mask%map(:,1)))
           mean = sum(hgm%map(:, imap), mask = mask%map(:,1) .ne. 0.)/weight
           rms = sqrt(sum((hgm%map(:, imap)-mean)**2, mask = mask%map(:,1) .ne. 0.)/weight)
           mapmin = minval(hgm%map(:,imap), mask = mask%map(:,1) .ne. 0.)
           mapmax = maxval(hgm%map(:,imap), mask = mask%map(:,1) .ne. 0.)
        else
           weight = hgm%npix
           mean = sum(hgm%map(:, imap))/weight
           rms = sqrt(sum((hgm%map(:, imap)-mean)**2)/weight)
           mapmin = minval(hgm%map(:,imap))
           mapmax =  maxval(hgm%map(:,imap))
        endif

        write(*, "(A)") "--- map #"//COOP_STR_OF(imap)//" --- "
        write(*, "(A, I5)") "spin = ", hgm%spin(imap)
        write(*, "(A)") "    #pix>0 fsky= "//COOP_STR_OF(count(hgm%map(:, imap).gt. 0.)/dble(hgm%npix))
        write(*, "(A)") "    #pix<0 fsky= "//COOP_STR_OF(count(hgm%map(:, imap).lt. 0.)/dble(hgm%npix))
        write(*, "(A)") "    min = "//COOP_STR_OF(mapmin)
        write(*, "(A)") "    max = "//COOP_STR_OF(mapmax)
        write(*, "(A)") "    mean = "//COOP_STR_OF(mean)
        write(*, "(A)") "    rms = "//COOP_STR_OF(rms)
     enddo
     write(*,*) "Print the HEADER (Y/N)?"
     read(*,*) yn
     if(yn .eq. "y" .or. yn .eq. "Y") &
          call hgm%header%print()
     call hgm%free()
  else
     write(*,"(A)") trim(fin)//" does not exist"
  endif
#else
  stop "you need to install healpix"
#endif
end program map
