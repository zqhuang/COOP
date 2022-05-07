program udg
  use coop_wrapper_utils
  use coop_healpix_mod
#ifdef HAS_HEALPIX  
  use udgrade_nr
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
#endif  
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX  
  character(LEN=1024)::f_map_in, f_map_out, f_mask, line
  integer nside_in, nside_out, npix_in, npix_out
  logical::ave, has_mask
  type(coop_healpix_maps)::map_in, map_out, mask
  if(iargc().lt.3)then
     write(*,*) "Syntax:"
     write(*,*) "./UDG -inp MAP_IN -out MAP_OUT -nside NSIDE_OUT -mask MASK -ave DO_AVERAGE[T/F]"
     write(*,*) "If -ave set to be false, will do mask udgrade (>0.5 ~ 1; <=0.5 ~ 0)."
     stop
  endif
  call coop_get_command_line_argument(key = "inp", arg  = f_map_in)
  call coop_get_command_line_argument(key = "out", arg = f_map_out)
  call coop_get_command_line_argument(key = "nside", arg = nside_out)
  call coop_get_command_line_argument(key = "ave", arg = ave, default = .true.)
  call coop_get_command_line_argument(key = "mask", arg = f_mask, default = "NONE")
  has_mask = (trim(f_mask) .ne. "NONE")
  call map_in%read(f_map_in)
  call map_out%init(nside = nside_out, nmaps = map_in%nmaps)
  call coop_healpix_maps_copy_genre(map_in, map_out)  
  if(has_mask)then
     call mask%read(f_mask)
     call map_in%apply_mask(mask)
     call coop_healpix_maps_ave_udgrade(from = map_in, to = map_out, mask = mask)     
  else
     call coop_healpix_maps_ave_udgrade(map_in, map_out)
  endif

  if(.not. ave .and. map_in%nmaps.eq.1 .and. all(map_in%map .eq. 0. .or. map_in%map .eq. 1.))then
     where(map_out%map .ge. 0.5)
        map_out%map = 1.
     elsewhere
        map_out%map = 0.
     end where
  endif
  !unit conversion
  if(trim(map_in%units(1)) .eq. "1")then
     if(coop_strpos(f_map_in, "_K_") .gt. 0 )then !!WMAP_K unit mK
        map_out%map = map_out%map * 1.e3
     elseif(coop_strpos(f_map_in, "HFI") .gt. 0)then !!Planck map unit K
        map_out%map = map_out%map * 1.e6
     endif
     call map_out%set_units("uK")
  endif
  if(map_in%nmaps .eq. 3 .or. map_in%nmaps .eq. 1)then     
     call map_out%write(trim(f_map_out))
  else if(map_in%nmaps .eq. 4)then
     call map_out%set_field(1, "T")
     call map_out%set_field(3, "Q")
     call map_out%set_field(4, "U")
     map_out%spin = (/ 0, 0, 2, 2 /)
     call map_out%write(trim(f_map_out), index_list = (/ 1, 3, 4 /) )
  endif
#else
  stop "you need to install healpix"
#endif  
end program udg
