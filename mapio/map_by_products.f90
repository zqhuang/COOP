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
  COOP_STRING::imap, qumap, tqtut, emap, bmap, zeta, zetaqzuz, outi, outqu, imask, polmask, str_fwhm_in, str_fwhm_out, prefix
  type(coop_healpix_maps) hgm, hgimask, hgpolmask
  COOP_INT fwhm_in_arcmin, fwhm_out_arcmin, lmax
  COOP_REAL fwhm_in, fwhm_out
  !!make code faster 
  coop_healpix_alm_check_done = .true.
  coop_healpix_want_cls = .false.
  if(iargc() < 7)then
     write(*,*) "Syntax:"
     write(*,*) "ByProd output_prefix imap polmap imask polmask fwhm_in_acrmin fwhm_out_arcmin"
     stop 
  endif
  prefix = trim(coop_inputArgs(1))
  imap = trim(coop_inputArgs(2))
  qumap = trim(coop_inputArgs(3))
  imask = trim(coop_inputArgs(4))
  polmask = trim(coop_inputArgs(5))
  str_fwhm_in = trim(coop_inputArgs(6))  
  str_fwhm_out = trim(coop_inputArgs(7))
  read(str_fwhm_in, *) fwhm_in_arcmin
  read(str_fwhm_out, *) fwhm_out_arcmin
  fwhm_in = fwhm_in_arcmin*coop_SI_arcmin
  fwhm_out = fwhm_out_arcmin*coop_SI_arcmin
  lmax = min(ceiling(3.d0/(fwhm_out*coop_sigma_by_fwhm)), 2500)
  if(.not. (coop_file_exists(trim(imap)) .and. coop_file_exists(trim(qumap))))then
     print*, "imap or qumap does not exist"
     stop
  endif
  tqtut = trim(prefix)//"_TQTUT_fwhm"//trim(str_fwhm_out)//".fits"
  zeta = trim(prefix)//"_zeta_fwhm"//trim(str_fwhm_out)//".fits"
  zetaqzuz = trim(prefix)//"_zetaqzuz_fwhm"//trim(str_fwhm_out)//".fits"  
  outi = trim(prefix)//"_I_fwhm"//trim(str_fwhm_out)//".fits"
  outqu = trim(prefix)//"_QU_fwhm"//trim(str_fwhm_out)//".fits"
  emap =trim(prefix)//"_E_fwhm"//trim(str_fwhm_out)//".fits"
  bmap = trim(prefix)//"_B_fwhm"//trim(str_fwhm_out)//".fits"
  if(.not. coop_file_exists(trim(outi)) .or. .not. coop_file_exists(trim(tqtut)) .or. .not. coop_file_exists(trim(zeta)) .or. .not. coop_file_exists(trim(zetaqzuz)) )then
     call hgm%read(trim(imap), nmaps_wanted = 3, spin=(/ 0, 2, 2 /), nmaps_to_read = 1 )
     if(trim(imask).ne."NONE")then
        call hgimask%read(trim(imask), nmaps_wanted = 1, spin = (/ 0 /))
        call hgm%regularize_in_mask(hgimask, 1)
     endif
     if(fwhm_out.gt.fwhm_in) &
          call hgm%smooth(fwhm = sqrt(fwhm_out**2-fwhm_in**2), index_list = (/ 1 /), l_upper = lmax)
     if(.not. coop_file_exists(trim(outi)))call hgm%write(trim(outi), index_list = (/ 1 /))
     if(.not. coop_file_exists(trim(tqtut)))then
        call hgm%iqu2TQTUT()
        call hgm%write(trim(tqtut) )
     endif
     if(.not. coop_file_exists(trim(zeta)) .or. .not. coop_file_exists(trim(zetaqzuz)) )then
        call hgm%t2zeta(fwhm_arcmin = 7.d0)  !!assumes default planck noise level
        if(.not. coop_file_exists(trim(zeta)))call hgm%write(trim(zeta), index_list = (/ 1 /) )
        if( .not. coop_file_exists(trim(zetaqzuz)))then
           call hgm%iqu2TQTUT()
           call hgm%write(trim(zetaqzuz))
        endif
     endif
  endif

  if(.not. coop_file_exists(trim(outqu)) .or. .not. coop_file_exists(trim(emap)) .or. .not. coop_file_exists(trim(bmap)))then
     call hgm%read(trim(qumap), nmaps_wanted = 2, spin = (/ 2, 2 /))
     if(trim(polmask).ne."NONE")then
        call hgpolmask%read(trim(polmask),nmaps_wanted = 1, spin = (/ 0 /) )
        call hgm%regularize_in_mask(hgpolmask, 1)
        call hgm%regularize_in_mask(hgpolmask, 2)        
     endif
     if(fwhm_out .gt. fwhm_in) &
          call hgm%smooth(fwhm = sqrt(fwhm_out**2-fwhm_in**2), l_upper = lmax)
     if(.not. coop_file_exists(trim(outqu)))call hgm%write(trim(outqu))
     if( .not. coop_file_exists(trim(emap)) .or.  .not. coop_file_exists(trim(bmap)))then
        call hgm%qu2EB()
        call hgm%write(trim(emap), index_list = (/ 1 /) )
        call hgm%write(trim(bmap), index_list = (/ 2 /) )
     endif
  endif
#else
  stop "You need to install healpix"
#endif  

  
end program map
