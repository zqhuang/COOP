program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  COOP_INT::lmin, lmax, lx_cut, ly_cut
  COOP_REAL fwhm_arcmin
  COOP_STRING::map, output, beamfile, mask
  COOP_SHORT_STRING::field, unit, genre
  type(coop_flatsky_maps)::fm
  type(coop_file) fp
  logical want_qu
  COOP_INT::l, i
  COOP_REAL,dimension(:),allocatable::beam
  COOP_REAL::reg_limit
  if(iargc().lt.2)then
     write(*,*) "Syntax:"
     write(*,*) "./FSmooth -map MAP_FILE -out OUTPUT  [-mask MASK] [-reg REGULARIZE_LIMIT] [-beam BEAM_FILE] [-fwhm FWHM_IN_ARCMIN ] [-lmin LMIN] [-lmax LMAX] [-lxcut LXCUT] [-lycut LYCUT] [-wantqu WANT_QU] [-field FIELD] [-unit UNIT]"     
     stop
  endif
  call coop_get_command_line_argument(key = "map", arg = map)
  call coop_get_command_line_argument(key = "out", arg = output)
  call coop_get_command_line_argument(key = "unit", arg = unit, default="muK")
  call coop_get_command_line_argument(key = "field", arg = field, default="T")
  call coop_get_command_line_argument(key = "wantqu", arg = want_qu,default=.false.)
  call coop_get_command_line_argument(key = "reg", arg = reg_limit, default=0.d0)
  call coop_get_command_line_argument(key = "mask", arg = mask,default="")
  call coop_get_command_line_argument(key = "beam", arg = beamfile, default="")
  call coop_get_command_line_argument(key = "fwhm", arg = fwhm_arcmin, default=0.d0)
  call coop_get_command_line_argument(key = "lmin", arg = lmin, default=2)
  call coop_get_command_line_argument(key = "lmax", arg = lmax, default=10000)
  call coop_get_command_line_argument(key = "lxcut", arg = lx_cut, default=0)
  call coop_get_command_line_argument(key = "lycut", arg = ly_cut, default=0)
  allocate(beam(0:lmax))
  if(trim(beamfile).ne."" .and. trim(beamfile).ne."NULL")then
     call fp%open_skip_comments(beamfile)
     do l = 0, lmax
        read(fp%unit, *) i, beam(l)
        if(i.ne.l) stop "beam file error"
     enddo
     call fp%close()
  else
     beam = 1.d0
  endif
  if(want_qu)then
     select case(COOP_UPPER_STR(field))
     case("T", "I")
        genre = "TQTUT"
     case("ZETA", "Z")
        genre = "ZQZUZ"
     case default
        stop "Error: wantqu is only enabled for T or zeta map"
     end select
  endif
  if(trim(mask).ne."")then
     if(want_qu)then
        call fm%read_from_one(filename = map, mask = mask, nmaps = 3, genre="TQTUT")
     else
        call fm%read_from_one(filename = map, mask = mask, nmaps = 1, genre=field)
     endif
  else
     if(want_qu)then
        call fm%read_from_one(filename = map, nmaps = 3, genre="TQTUT")
     else
        call fm%read_from_one(filename = map, nmaps = 1, genre=field)
     endif
  endif
  call fm%map(1)%regularize(reg_limit)
  
  call fm%map(1)%smooth(fwhm = fwhm_arcmin*coop_SI_arcmin, highpass_l1 = max(lmin-20, 2), highpass_l2 = lmin + 20, lmax = lmax, lx_cut = lx_cut, ly_cut = ly_cut, beam = beam)
  fm%units = unit
  if(want_qu)call fm%map(1)%get_QTUT(qt = fm%map(2), ut = fm%map(3))
  call fm%write(output)
  write(*,*) "File is written to "//trim(output)

end program test
