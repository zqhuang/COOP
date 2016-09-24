program fsm
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  COOP_INT::lmin, lmax, lx_cut, ly_cut
  COOP_REAL fwhm_arcmin
  COOP_STRING::qmap, umap,  beamfile, mask, out_prefix
  COOP_SHORT_STRING::field, unit, genre
  type(coop_flatsky_maps)::fm,eb
  type(coop_file) fp
  COOP_INT::l, i
  COOP_REAL,dimension(:),allocatable::beam
  COOP_REAL::reg_limit
  logical::overwrite
  if(iargc().lt.2)then
     write(*,*) "Syntax:"
     write(*,*) "./FSmooth -qmap ... -umap  -out_prefix ... [-mask MASK] [-reg REGULARIZE_LIMIT] [-beam BEAM_FILE] [-fwhm FWHM_IN_ARCMIN ] [-lmin LMIN] [-lmax LMAX] [-lxcut LXCUT] [-lycut LYCUT] [-field FIELD] [-unit UNIT] [-overwrite F/T]"     
     stop
  endif
  call coop_get_command_line_argument(key = "out_prefix", arg = out_prefix)
  call coop_get_command_line_argument(key = 'overwrite', arg= overwrite, default = .true.)
  if(coop_file_exists(trim(out_prefix)//"_Q.fsm") .and. coop_file_exists(trim(out_prefix)//"_U.fsm") .and. coop_file_exists(trim(out_prefix)//"_E.fsm") .and. coop_file_exists(trim(out_prefix)//"_B.fsm") .and. .not. overwrite)then
     write(*, *) "the output files "//trim(out_prefix)//"_Q/U/E/B.fsm already exist"
  else
     call coop_get_command_line_argument(key = "qmap", arg = qmap)
     call coop_get_command_line_argument(key = "umap", arg = umap)
     call coop_get_command_line_argument(key = "unit", arg = unit, default="muK")
     call coop_get_command_line_argument(key = "field", arg = field, default="T")
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
     if(trim(mask).ne."")then
        call fm%read_from_one(filename = qmap, mask = mask, nmaps = 2, genre="QU")
        call eb%read_from_one(filename = qmap, mask = mask, nmaps = 2, genre="EB")
     else
        call fm%read_from_one(filename = qmap, nmaps = 2, genre="QU")
        call eb%read_from_one(filename = qmap, nmaps = 2, genre="EB")
     endif
     call fm%map(2)%open(umap)
     if(reg_limit .gt. 0.d0)then
        call fm%map(1)%regularize(reg_limit)
        call fm%map(2)%regularize(reg_limit)
     endif
     if(fm%has_mask .and. fm%mask_is_smooth)then
        write(*,* ) "warning: applying mask before smoothing and QU->EB"
        where(.not. fm%unmasked)
           fm%map(1)%image = fm%map(1)%image * sin(fm%mask%image/fm%mask_threshold*coop_pio2)**2
        end where
     endif
     call coop_fits_image_cea_smooth_qu(qmap = fm%map(1), umap=fm%map(2), emap = eb%map(1), bmap=eb%map(2), fwhm = fwhm_arcmin*coop_SI_arcmin, highpass_l1 = max(lmin-20, 2), highpass_l2 = lmin + 20, lmax = lmax, lx_cut = lx_cut, ly_cut = ly_cut, beam = beam)
     fm%map_changed = .true.
     eb%map_changed = .true.
     fm%units = unit
     eb%units = unit
     call fm%write(trim(out_prefix)//"_QU.fsm")
     call fm%write(trim(out_prefix)//"_Q.fsm", write_image = .false., index_list = (/ 1 /) )
     call fm%write(trim(out_prefix)//"_U.fsm",  write_image = .false., index_list = (/ 2 /) )
     call eb%write(trim(out_prefix)//"_EB.fsm")
     call eb%write(trim(out_prefix)//"_E.fsm", write_image = .false., index_list = (/ 1 /) )
     call eb%write(trim(out_prefix)//"_B.fsm",  write_image = .false., index_list = (/ 2 /) )
     write(*,*) "Files saved: "//trim(out_prefix)//"_Q/U/E/B.fsm"
  endif
end program fsm
