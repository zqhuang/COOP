program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  type(coop_fits_image_cea)::total_map, total_weights, this_map, this_weights
  COOP_STRING::params_file, map_file, weight_file
  logical::positive_weights = .true.
logical::analyze_maps = .false.
  type(coop_dictionary)::params
  COOP_INT:: num_maps, i
  COOP_REAL::coef, truncate, mean_weight
!  call coop_MPI_Init()
  call coop_get_Input(1, params_file)
  call coop_load_dictionary(params_file, params)
  call coop_dictionary_lookup(params, "num_maps", num_maps)
  if(num_maps .gt. 100 .or. num_maps.lt.1)then
     stop "num_maps must be between 1 and 100"
  else
     write(*,*) "coadding "//COOP_STR_OF(num_maps)//" maps"
  endif
  call coop_dictionary_lookup(params, "positive_weights", positive_weights, default_val = .true.)
  call coop_dictionary_lookup(params, "analyze_maps", analyze_maps, default_val = .false.)
  do i=1, num_maps
     call coop_dictionary_lookup(params, "map"//COOP_STR_OF(i), map_file)
     if(.not. coop_file_exists(map_file))then
        write(*,*) "map"//COOP_STR_OF(i)//":"//trim(map_file)//" is missing"
        stop
     endif
     call coop_dictionary_lookup(params, "weight"//COOP_STR_OF(i), weight_file)
     if(.not. coop_file_exists(weight_file))then
        write(*,*) "weight"//COOP_STR_OF(i)//":"//trim(weight_file)//" is missing"
        stop
     endif
     call coop_dictionary_lookup(params, "coef"//COOP_STR_OF(i), coef, default_val = 1.d0)
     if(i.eq.1)then
        call total_map%open(map_file)
        if(analyze_maps)call total_map%simple_stat()
        call total_weights%open(weight_file)
        if(analyze_maps)call total_weights%simple_stat()
        if(total_map%npix .ne. total_weights%npix)then
           write(*,*) "Error: map and weight must be the same size"
           stop
        endif
        if(positive_weights)then
           total_weights%image = abs(total_weights%image)*coef !!only positive values
        elseif(abs(coef-1.d0) .gt. 1.d-6)then
           total_weights%image = total_weights%image*coef
        endif
        total_map%image = total_map%image*total_weights%image
     else
        call this_map%open(map_file)
        if(analyze_maps)call this_map%simple_stat()
        call this_weights%open(weight_file)
        if(analyze_maps)call this_weights%simple_stat()
        if(this_map%npix .ne. total_map%npix .or. this_weights%npix .ne. total_weights%npix)then
           write(*,*) "maps with different sizes cannot be coadded"
           stop
        endif
        if(positive_weights)then
           this_weights%image = abs(this_weights%image)*coef !!only positive values
        elseif(abs(coef-1.d0) .gt. 1.d-6)then
           this_weights%image= this_weights%image*coef
        endif
        total_map%image = total_map%image + this_map%image*this_weights%image
        total_weights%image = total_weights%image + this_weights%image
     endif
  enddo
  call this_map%free()
  call this_weights%free()
  call coop_dictionary_lookup(params, "truncate_weight", truncate, default_val = 0.05d0)
  mean_weight = sum(total_weights%image)/total_weights%npix
  truncate = truncate*mean_weight
  write(*,*) "mean weight = ", mean_weight
  write(*,*) "truncating at weight <=", truncate
  where(total_weights%image .gt. truncate)
     total_map%image = total_map%image/total_weights%image
  elsewhere 
     total_map%image = 0.d0
     total_weights%image = 0.d0
  end where
  write(*,*) "# pixels = ",total_map%npix
  write(*,*) "fraction of truncated pixels:", count(total_weights%image .eq. 0.d0)/dble(total_weights%npix)
  call coop_dictionary_lookup(params, "output_map", map_file)
  call total_map%write(map_file)
  call coop_dictionary_lookup(params, "output_weight", weight_file)
  call total_weights%write(weight_file)
end program test
