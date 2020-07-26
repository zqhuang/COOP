program Stacking_Maps
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX
  logical::remove_mono = .true.
  logical::randrot, want_pdf
  COOP_STRING::mask_file, peaks_file, map_file, output
  COOP_REAL::r_degree, dr
  type(coop_stacking_options)::sto
  type(coop_healpix_filament)::filament
  type(coop_healpix_maps)::hgm, mask
  COOP_INT i, m, ny, n
  COOP_REAL::zmin1, zmax1, zmin2, zmax2
  COOP_REAL::tmax
  call random_seed()
  if(iargc() .lt. 6) then
     write(*,*) "./StackFil -peaks peaks_file -map map_file -out output_file"
     stop
  endif
  call coop_get_command_line_argument(key = 'peaks', arg = peaks_file)
  call coop_get_command_line_argument(key = 'map', arg = map_file)
  call coop_get_command_line_argument(key = 'out', arg = output)
  call coop_get_command_line_argument(key = 'width', arg = r_degree, default = 7.d0)
  call coop_get_command_line_argument(key = 'resx', arg = n, default = 100)
  call coop_get_command_line_argument(key = 'resy', arg = ny, default = n)  
  coop_healpix_patch_default_figure_height = 3.5*ny/n  + 0.5
  dr = r_degree*coop_SI_degree/n
  call sto%import(peaks_file)
  call hgm%read(map_file)
  call filament%init(genre = "T", n = n, dr = dr, ny = ny)
  
  print*, "stacking on "//COOP_STR_OF(sto%peak_pix%n)//" peaks"  
  call hgm%stack_filaments_on_peaks(sto, filament)
  print*, "done"
  filament%caption = "stacked on "//COOP_STR_OF(sto%peak_pix%n)//" "//trim(sto%caption)
  call filament%plot(1, trim(adjustl(output)))
  call system("../utils/fasy.sh "//trim(adjustl(output)))
#endif
  
end program Stacking_Maps
