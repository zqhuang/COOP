program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  type(coop_dictionary)::params
  type(coop_fits_image_cea)::map, map2, mask, mask2
  logical::print_header = .false.
  COOP_STRING::fname, fname2, fmask, fmask2, cl_file
  COOP_REAL::mean, rms, mean2, rms2, threshold, cut, reg_limit

  if(iargc().lt.1)then
     write(*,*) "Syntax:"
     write(*,*) "./FStat -map MAPFILE [-mask MASKFILE] [-map2 SECOND_MAPFILE] [-mask2 SECOND_MASKFILE] [-cut CUT(0.1)] [-reg_limit REG_LIMIT(0)] [-cl_file CL_FILE]"
     stop
  endif
  if(iargc().eq.1)then
     fname = trim(coop_InputArgs(1))
     fname2 = ''
     fmask = ''
     cl_file = ''
     fmask2 = ''
     reg_limit = 0.d0
     cut = 0.1d0
  else
     call coop_get_command_line_argument(key = 'map', arg = fname)
     call coop_get_command_line_argument(key = 'print_header', arg = print_header, default = .false.)
     call coop_get_command_line_argument(key = 'map2', arg = fname2, default="")
     call coop_get_command_line_argument(key = 'cl_file', arg = cl_file, default="")
     call coop_get_command_line_argument(key = 'mask', arg = fmask, default="")
     call coop_get_command_line_argument(key = 'mask2', arg = fmask2, default="")
     call coop_get_command_line_argument(key = 'cut', arg = cut, default = 0.1d0)
     call coop_get_command_line_argument(key = 'reg_limit', arg = reg_limit, default = 0.d0)
  endif
  call map%open(fname)
  call map%regularize(reg_limit)
  if(trim(fmask).ne."")then
     call mask%open(fmask)
     if(map%npix .ne. mask%npix) stop "map and mask are of different sizes"
     threshold = maxval(mask%image)*cut
     where(mask%image .lt. threshold)
        mask%image = 0.
        map%image = 0.
     elsewhere
        mask%image = 1.
     end where
  endif
  if(trim(fname2).ne."")then
     call map2%open(fname2)
     call map2%regularize(reg_limit)
     if(map2%npix .ne. map%npix) stop "two maps are of different sizes"
     if(trim(fmask2).ne."")then
        call mask2%open(fmask2)
        if(map2%npix .ne. mask2%npix) stop "map and mask are of different sizes"
        threshold = maxval(mask2%image)*cut
        where(mask2%image .lt. threshold)
           mask2%image = 0.
           map2%image = 0.
        elsewhere
           mask2%image = 1.
        end where
        if(trim(fmask).ne."")then
           mask%image = mask%image*mask2%image
        else
           mask = mask2
        endif
     endif
  endif
  if(trim(fmask).ne."" .or. trim(fmask2).ne."")then
     call map%simple_stat(mean=mean, rms=rms, mask=mask,clsfile=cl_file)
  else
     call map%simple_stat(mean=mean, rms=rms,clsfile=cl_file)
  endif
  if(trim(fname2).ne."")then
     if(trim(fmask).ne."" .or. trim(fmask2).ne."")then
        call map2%simple_stat(mean=mean2, rms=rms2, mask = mask)
        write(*,*) "correlation between two maps:", sum((map%image-mean)*(map2%image-mean2), mask = mask%image .gt. 0.d0)/sqrt(sum((map%image-mean)**2, mask = mask%image .gt. 0.d0)*sum((map2%image-mean2)**2, mask = mask%image .gt. 0.d0))
     else
        call map2%simple_stat(mean=mean2, rms=rms2)
        write(*,*) "correlation between two maps:", sum((map%image-mean)*(map2%image-mean2))/sqrt(sum((map%image-mean)**2)*sum((map2%image-mean2)**2))


     endif
  endif
  if(print_header)then
     write(*,*) "============= header of "//trim(fname)//" =========="
     call map%header%print()
     write(*,*) "===================================================="
     if(trim(fname2) .ne. '')then
        write(*,*) "============= header of "//trim(fname2)//" =========="
        call map2%header%print()
        write(*,*) "===================================================="
     endif
  endif
  call map%free()
  call map2%free()
end program test
