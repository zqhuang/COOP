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
  
  COOP_STRING::imap, qumap, iqtut, emap, bmap
  type(coop_healpix_maps) hgm
  imap = trim(coop_inputArgs(1))
  qumap = trim(coop_inputArgs(2))
  if(trim(imap).eq."" .or. trim(qumap).eq."")then
     print*, "ByProd imap qumap"
     stop 
  endif
  if(.not. (coop_file_exists(trim(imap)) .and. coop_file_exists(trim(qumap))))then
     print*, "imap or qumap does not exist"
     stop
  endif
  iqtut = trim(coop_file_add_postfix(trim(imap), "_converted_to_TQTUT"))

  if(.not. coop_file_exists(trim(iqtut)))then
     call hgm%read(trim(imap), nmaps_wanted = 3, spin=(/ 0, 2, 2 /), nmaps_to_read = 1 )
     call hgm%iqu2TQTUT()
     call hgm%write(trim(iqtut))
  endif
  emap = trim(coop_file_add_postfix(qumap, "_converted_to_EB_E"))
  bmap = trim(coop_file_add_postfix(qumap, "_converted_to_EB_B"))
  if(.not.coop_file_exists(trim(emap)) .or. .not. coop_file_exists(trim(bmap)))then
     call hgm%read(trim(qumap), nmaps_wanted = 2, spin = (/ 2, 2 /))
     call hgm%qu2EB()
     if(.not.coop_file_exists(trim(emap)))call hgm%write(trim(emap), index_list = (/ 1 /) )
     if(.not. coop_file_exists(trim(bmap)))call hgm%write(trim(bmap), index_list = (/ 2 /) )
  endif
end program map
