program flatcoadd
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
  logical::do_filtering = .false.
  logical::has_weights = .true.
  type(coop_dictionary)::params
  COOP_INT:: num_maps, i, lmin, lmax
  COOP_REAL::coef, truncate, mean_weight, fwhm, reg_limit
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
  call coop_dictionary_lookup(params, "do_filtering", do_filtering, .false.)
  if(do_filtering)then
     write(*,*) "Doing filtering before coadding."
  endif
  call coop_dictionary_lookup(params, "highpass_lmin", lmin, default_val = 100)
  call coop_dictionary_lookup(params, "lowpass_lmax", lmax, default_val = 4000)
  call coop_dictionary_lookup(params, "fwhm_arcmin", fwhm, default_val = 0.d0)
  call coop_dictionary_lookup(params, "reg_limit", reg_limit, default_val = 0.d0)
  fwhm =fwhm*coop_SI_arcmin

  has_weights = .true.
  do i=1, num_maps
     call coop_dictionary_lookup(params, "map"//COOP_STR_OF(i), map_file)
     if(.not. coop_file_exists(map_file))then
        write(*,*) "map"//COOP_STR_OF(i)//":"//trim(map_file)//" is missing"
        stop
     endif
     if(has_weights)then
        call coop_dictionary_lookup(params, "weight"//COOP_STR_OF(i), weight_file)
        if(trim(weight_file).eq. "")then
           if(i.eq.1)then
              has_weights = .false.
              write(*,*) "Warning: weight"//COOP_STR_OF(i)//" is not set; using uniform weights."
           else
              write(*,*) "Error: weight"//COOP_STR_OF(i)//" is not set."
              stop
           endif
        else
           if(.not. coop_file_exists(weight_file))then
              write(*,*) "weight"//COOP_STR_OF(i)//":"//trim(weight_file)//" is missing"
              stop
           endif
        endif
     endif
     call coop_dictionary_lookup(params, "coef"//COOP_STR_OF(i), coef, default_val = 1.d0)
     if(i.eq.1)then
        call total_map%open(map_file)
        call total_map%regularize(reg_limit)
        if(do_filtering)call total_map%smooth(fwhm = fwhm, highpass_l1 = lmin - 10, highpass_l2 = lmin + 10, lowpass_l1 = lmax - 100, lowpass_l2 = lmax+100)
        if(analyze_maps)call total_map%simple_stat()
        if(has_weights)then
           call total_weights%open(weight_file)
           if(analyze_maps)call total_weights%simple_stat()
           if(positive_weights)then
              if(any(total_weights%image .lt. 0.d0)) stop "found negative weights"
           endif
           if(total_map%npix .ne. total_weights%npix)then
              write(*,*) "maps/mask with different sizes cannot be coadded"
              stop
           endif
           if(abs(coef-1.d0) .gt. 1.d-6)then
              total_map%image = total_map%image*total_weights%image*coef
           else
              total_map%image = total_map%image*total_weights%image
           endif
        elseif(abs(coef-1.d0) .gt. 1.d-6)then
           total_map%image = total_map%image * coef
        endif
     else
        call this_map%open(map_file)
        call this_map%regularize(reg_limit)
        if(do_filtering)call this_map%smooth(fwhm = fwhm, highpass_l1 = lmin - 10, highpass_l2 = lmin + 10, lowpass_l1 = lmax - 100, lowpass_l2 = lmax+100)
        if(analyze_maps)call this_map%simple_stat()
        if(has_weights)then
           call this_weights%open(weight_file)
           if(analyze_maps)call this_weights%simple_stat()
           if(positive_weights)then
              if(any(this_weights%image .lt. 0.d0)) stop "found negative weights"
           endif
           if(this_map%npix .ne. total_map%npix .or. this_weights%npix .ne. total_weights%npix)then
              write(*,*) "maps with different sizes cannot be coadded"
              stop
           endif
           if(abs(coef-1.d0) .gt. 1.d-6)then
              total_map%image = total_map%image + this_map%image*this_weights%image*coef
           else
              total_map%image = total_map%image + this_map%image*this_weights%image
           endif
           total_weights%image = total_weights%image + this_weights%image
        else
           if(abs(coef-1.d0) .gt. 1.d-6)then
              total_map%image = total_map%image + this_map%image * coef
           else
              total_map%image = total_map%image + this_map%image
           endif
        endif
     endif
  enddo
  call this_map%free()
  call this_weights%free()
  call coop_dictionary_lookup(params, "output_map", map_file)
  call coop_dictionary_lookup(params, "output_weight", weight_file)
  if(trim(map_file).eq."") stop "output_map is not set"
  if(has_weights)then
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
     write(*,*) "fraction of truncated pixels:", count(total_weights%image .eq. 0.d0)/dble(total_weights%npix)
     if(trim(weight_file).ne."") call total_weights%write(weight_file)
     call total_map%write(map_file)
  else
     total_map%image = total_map%image/num_maps
     call total_map%write(map_file)
     if(trim(weight_file).ne."")then
        total_map%image = 1.d0
        call total_map%write(weight_file)
     endif
  endif
  write(*,*) "The coadded map has been written to "//trim(map_file)
  call total_map%free()
  call total_weights%free()
end program flatcoadd


