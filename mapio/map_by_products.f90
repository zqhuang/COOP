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
  
  COOP_STRING::imap, qumap, iqtut, emap, bmap, zeta, imask, polmask
  type(coop_healpix_maps) hgm, hgimask, hgpolmask
  COOP_REAL fwhm_arcmin
  !!make code faster 
  coop_healpix_alm_check_done = .true.
  coop_healpix_want_cls = .false.

  imap = trim(coop_inputArgs(1))
  qumap = trim(coop_inputArgs(2))
  zeta = trim(coop_inputArgs(3))
  imask = trim(coop_inputArgs(4))
  polmask = trim(coop_inputArgs(5))

  if(trim(zeta).eq. "")then
     fwhm_arcmin = 10.d0  !!default planck nside 1024
  else
     read(zeta, *) fwhm_arcmin
     if(fwhm_arcmin .lt. 0.d0 .or. fwhm_arcmin .gt. 1000.d0)then
        print*, "ByProd imap qumap fwhm_arcmin"
        stop "fwhm_arcmin not in the right range (0 - 1000)"
     endif
  endif
  if(trim(imap).eq."" .or. trim(qumap).eq."")then
     print*, "ByProd imap qumap fwhm_arcmin"
     stop 
  endif
  if(.not. (coop_file_exists(trim(imap)) .and. coop_file_exists(trim(qumap))))then
     print*, "imap or qumap does not exist"
     stop
  endif
  iqtut = trim(coop_file_add_postfix(trim(imap), "_converted_to_TQTUT"))
  zeta = trim(coop_file_add_postfix(trim(imap), "_converted_to_zeta"))

  if(.not. coop_file_exists(trim(iqtut)) .or. .not. coop_file_exists(trim(zeta)) )then
     call hgm%read(trim(imap), nmaps_wanted = 3, spin=(/ 0, 2, 2 /), nmaps_to_read = 1 )
     if(trim(imask).ne."")then
        call hgimask%read(trim(imask), nmaps_wanted = 1)
        hgm%map(:,1) = hgm%map(:,1)*hgimask%map(:, 1)
     endif
     if(.not. coop_file_exists(trim(iqtut)))then
        call hgm%iqu2TQTUT()
        call hgm%write(trim(iqtut))
     endif
     if(.not. coop_file_exists(trim(zeta)))then
        call hgm%t2zeta( fwhm_arcmin)
        call hgm%write(trim(zeta), index_list = (/ 1 /) )
     endif
  endif



  emap = trim(coop_file_add_postfix(qumap, "_converted_to_EB_E"))
  bmap = trim(coop_file_add_postfix(qumap, "_converted_to_EB_B"))
  if( .not. coop_file_exists(trim(emap)) .or. .not. coop_file_exists(trim(bmap)))then
     call hgm%read(trim(qumap), nmaps_wanted = 2, spin = (/ 2, 2 /))
     if(trim(polmask).ne."")then
        call hgpolmask%read(trim(polmask),nmaps_wanted = 1)
        hgm%map(:, 1) = hgm%map(:, 1)*hgpolmask%map(:, 1)
        hgm%map(:, 2) = hgm%map(:, 2)*hgpolmask%map(:, 1)
     endif
     call hgm%qu2EB()
     if(.not.coop_file_exists(trim(emap)))call hgm%write(trim(emap), index_list = (/ 1 /) )
     if(.not. coop_file_exists(trim(bmap)))call hgm%write(trim(bmap), index_list = (/ 2 /) )
  endif

  

  
end program map
