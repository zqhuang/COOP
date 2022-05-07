program flatcoadd
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  type(coop_fits_image_cea)::total_map, total_weights, this_map, this_weights, map_neg, weight_neg
  COOP_STRING::params_file, map_file, weight_file, beam_file
  logical::positive_weights = .true.
  logical::analyze_maps = .false.
  logical::do_filtering = .false.
  logical::has_weights = .true.
  type(coop_dictionary)::params
  logical::do_diff = .false.
  COOP_INT:: num_maps, i, lmin, lmax, il, l
  COOP_REAL::coef, truncate, mean_weight, fwhm, reg_limit, truncate_weight
  COOP_REAL,dimension(:),allocatable::beam
  type(coop_file)::fp
  call coop_get_Input(1, params_file)
  call coop_load_dictionary(params_file, params)
  call coop_dictionary_lookup(params, "num_maps", num_maps)
  call coop_dictionary_lookup(params, "do_diff", do_diff, .false.)
  if(num_maps .gt. 100 .or. num_maps.lt.1)then
     stop "num_maps must be between 1 and 100"
  else
     write(*,*) "coadding "//COOP_STR_OF(num_maps)//" maps"
  endif
  call coop_dictionary_lookup(params, "positive_weights", positive_weights, default_val = .true.)
  call coop_dictionary_lookup(params, "analyze_maps", analyze_maps, default_val = .false.)
  call coop_dictionary_lookup(params, "reg_limit", reg_limit, default_val = 0.d0)
  call coop_dictionary_lookup(params, "do_filtering", do_filtering, .false.)
  if(do_filtering)then
     write(*,*) "Doing filtering before coadding."
     call coop_dictionary_lookup(params, "highpass_lmin", lmin, default_val = 100)
     write(*,*) "lmin = ", lmin
     call coop_dictionary_lookup(params, "lowpass_lmax", lmax, default_val = 4000)
     write(*,*) "lmax = ", lmax
     call coop_dictionary_lookup(params, "fwhm_arcmin", fwhm, default_val = 0.d0)
     write(*,*) "FWHM = ", nint(fwhm), " arcmin"
     fwhm =fwhm*coop_SI_arcmin
     allocate(beam(0:lmax))
  endif

  has_weights = .true.
  do i=1, num_maps
     write(*,*) "==========================================================="
     call coop_dictionary_lookup(params, "map"//COOP_STR_OF(i), map_file)
     write(*,*)  "Map file: "//trim(map_file)
     if(.not. coop_file_exists(map_file))then
        write(*,*) "map"//COOP_STR_OF(i)//":"//trim(map_file)//" is missing"
        stop
     endif
     if(has_weights)then
        call coop_dictionary_lookup(params, "weight"//COOP_STR_OF(i), weight_file)
        write(*,*)  "Weight file: "//trim(weight_file)
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
     call this_map%open(map_file)
     call this_map%regularize(reg_limit)
     if(do_filtering)then
        call coop_dictionary_lookup(params, "beam"//COOP_STR_OF(i), beam_file, "")
        if(trim(beam_file).ne."")then
           call fp%open_skip_comments(beam_file)
           do l = 0, lmax
              read(fp%unit, *) il, beam(l)
              if(il.ne.l) stop "beam file error"
           enddo
           call fp%close()   
           write(*,*) "debeaming using file "//trim(beam_file)
        else
           beam = 1.d0
        endif
        call this_map%smooth(fwhm = fwhm, highpass_l1 = lmin - 10, highpass_l2 = lmin + 10, lmax = lmax, beam = beam)
     endif
     if(analyze_maps)call this_map%simple_stat()
     if(i.eq.1)then
        total_map = this_map
        total_map%image = 0.
        if(has_weights)total_weights = total_map
        if(do_diff)then
           map_neg = total_map
           if(has_weights)weight_neg = total_weights
        endif
     endif
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
        if(do_diff .and. coef .lt. 0.)then
           map_neg%image = map_neg%image + this_map%image*this_weights%image*coef
           weight_neg%image = weight_neg%image + this_weights%image
        else
           total_map%image = total_map%image + this_map%image*this_weights%image*coef
           total_weights%image = total_weights%image + this_weights%image
        endif
     else
        total_map%image = total_map%image + this_map%image * coef
     endif
  enddo
  call this_map%free()
  call this_weights%free()
  call coop_dictionary_lookup(params, "output_map", map_file)
  call coop_dictionary_lookup(params, "output_weight", weight_file)
  if(trim(map_file).eq."") stop "output_map is not set"
  if(has_weights)then
     call coop_dictionary_lookup(params, "truncate_weight", truncate_weight, default_val = 0.05d0)
     mean_weight = sum(total_weights%image)/total_weights%npix
     truncate = truncate_weight * mean_weight
     write(*,*) "mean weight = ", mean_weight
     write(*,*) "truncating at weight <=", truncate
     where(total_weights%image .gt. truncate)
        total_map%image = total_map%image/total_weights%image
     elsewhere 
        total_map%image = 0.d0
        total_weights%image = 0.d0
     end where
     write(*,*) "fraction of truncated pixels:", count(total_weights%image .eq. 0.d0)/dble(total_weights%npix)
     if(do_diff)then
        mean_weight = sum(weight_neg%image)/weight_neg%npix
        truncate = truncate_weight*mean_weight
        where(weight_neg%image .gt. truncate)
           map_neg%image = map_neg%image/weight_neg%image
        elsewhere 
           map_neg%image = 0.d0
           weight_neg%image = 0.d0
        end where
        total_map%image = (total_map%image + map_neg%image)/2.d0
     endif
     call total_map%write(map_file)
     if(trim(weight_file).ne."")then
        if(do_diff)then
           total_weights%image = total_weights%image * weight_neg%image
        endif
        call total_weights%write(weight_file)
     endif
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
  call map_neg%free()
  call weight_neg%free()
end program flatcoadd


