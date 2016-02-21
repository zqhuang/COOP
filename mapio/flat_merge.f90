program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_flatsky_maps)::map1, map2, map3
  COOP_STRING::f1, f2, f3, fout
  logical mask_changed
  if(iargc().lt.4) stop "Syntax: ./FMerge -map1 MAP1 -map2 MAP2 [-map3 MAP3] -out OUTPUT"
  call coop_get_command_line_argument(key = "map1", arg = f1)
  call coop_get_command_line_argument(key = "map2", arg = f2)
  call coop_get_command_line_argument(key = "map3", arg = f3, default="")
  call coop_get_command_line_argument(key = "out", arg = fout)
  if(coop_file_exists(fout))then
     write(*,*) "the output file "//trim(fout)//" already exits"
  else
     call map1%read(f1)
     call map2%read(f2)
     call map1%merge(map2)
     if(trim(f3).ne."")then
        call map3%read(f3)
        call map1%merge(map3)
     endif
     call map1%write(fout, write_image = .false., write_mask = map1%mask_changed)
  endif
end program test
