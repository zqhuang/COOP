program Exp_spots
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX

  COOP_REAL::l_deg_axis = 220.d0
  COOP_REAL::b_deg_axis = -17.d0

  COOP_STRING::cc_method = "smica"
  COOP_STRING::output 
  COOP_REAL::threshold = 0.5d0
  
  COOP_INT,parameter::n_sim = 30
  COOP_STRING::peak_name = "$T$"
  COOP_STRING::orient_name = "NULL"
  COOP_STRING::stack_field_name = "T"
  COOP_INT,parameter::n = 36
  COOP_REAL,parameter::r_degree  = 2.d0
  COOP_REAL,parameter::dr = 2.d0*sin(r_degree*coop_SI_degree/2.d0)/n

  
  COOP_REAL::r(0:n), pfr(0:n), pfr(0:n)
  COOP_UNKNOWN_STRING,parameter::mapdir = "massffp8/"
  COOP_UNKNOWN_STRING,parameter::postfix = "_020a_0512.fits"
  COOP_STRING::imap_file = "planck14/dx11_v2_"//trim(cc_method)//"_int_cmb"//postfix
  COOP_STRING::polmap_file = "planck14/dx11_v2_"//trim(cc_method)//"_pol_case1_cmb_hp_20_40"//postfix  
  COOP_STRING::imask_file = "planck14/dx11_v2_common_int_mask"//postfix 
  COOP_STRING::polmask_file = "planck14/dx11_v2_common_pol_mask"//postfix
  COOP_STRING::line
  type(coop_stacking_options)::sto_max, sto_min, stn_max, sts_max, stn_min, sts_min
  type(coop_healpix_patch)::north_max, north_min, south_max, south_min
  type(coop_healpix_maps)::imap, imask, polmask, inoise, polnoise, polmap

  logical::iloaded = .false.
  logical::polloaded  = .false.
  type(coop_file)::fp
  type(coop_asy)::fig
  COOP_INT i
  COOP_INT ind, ind_done


  if(iargc().gt.3)then
     
  endif
  
  output = "ha/"//trim(cc_method)//"_hem_2deg_nu0pt5.dat"
  
  if(trim(orient_name).eq."NULL")then
     call sto_max%init(.true., peak_name, orient_name, nmaps = 1)
     call sto_min%init(.false., peak_name, orient_name, nmaps = 1)     
  else
     call sto_max%init(.true., peak_name, orient_name, nmaps = 3)          
     call sto_min%init(.false., peak_name, orient_name, nmaps = 3)     
  endif
  sto_max%I_lower_nu = threshold
  sto_min%I_upper_nu = -threshold
  sto_max%threshold_option = 4
  sto_min%threshold_option = 4
  sto_max%nested = .true.
  sto_min%nested = .true.  

  do i=0, n
     r(i) = i*dr
  enddo
    
  call north_max%init(trim(stack_field_name), n, dr)
  north_min = north_max
  south_max = north_max
  south_min = north_max
  
  ind_done = -1  
  call imask%read(imask_file, nmaps_wanted = 1, spin = (/ 0 /) )
  if(coop_file_exists(output))then
     call fp%open(output, "ru")
     do i = 0, n_sim
        read(fp%unit, ERR = 100, END = 100) ind, north_max%fr, north_min%fr, south_max%fr, south_min%fr
        ind_done = ind
     enddo
100  call fp%close()
  endif
  
  call fp%open(output, "u")
  call fig%open("fr_diff.txt")
  call fig%init(xlabel = "$\varpi$", ylabel = "$ T_0^{\rm north} - T_0^{\rm south}$")
  if(ind_done .ge. 0)then
     print*, "loaded "//COOP_STR_OF(ind_done+1)//" stacked maps"
  endif

  do ind = 0, n_sim
     if(ind.gt.ind_done)then
        print*, "stacking map#"//COOP_STR_OF(ind)                
        call load_imap(ind)
        call stack_imap()

        write(fp%unit) ind,  north_max%fr, north_min%fr, south_max%fr, south_min%fr
     else
        read(fp%unit) i, north_max%fr, north_min%fr, south_max%fr, south_min%fr
     endif
     pfr =  (north_max%fr(:, 0, 1) - north_min%fr(:,0,1))/2.d0 - (south_max%fr(:,0,1) - south_min%fr(:,0,1))/2.d0     
     if(ind.eq.0)then
        pfr0 = pfr
     elseif(ind.eq.1)then
        call fig%curve(r, pfr - pfr0, color = "blue", linetype= "dotted", linewidth = 1.,legend = "FFP8")        
     else
        call fig%curve(r, pfr - pfr0, color = "blue", linetype= "dotted", linewidth = 1.)
     endif
  enddo
  call coop_asy_legend(fig, "N", 2)
  call fig%close()
  call fp%close()
  

contains

  subroutine stack_imap()
    call imap%get_peaks(sto_max, mask = imask, restore = .false.)
    call imap%get_peaks(sto_min, mask = imask, restore = .false.)
    call coop_stacking_options_split_hemispheres(sto_max, stn_max, sts_max, l_deg_axis, b_deg_axis)
    call coop_stacking_options_split_hemispheres(sto_min, stn_min, sts_min, l_deg_axis, b_deg_axis)    

    call imap%stack_on_peaks(stn_max, north_max, imask)
    call imap%stack_on_peaks(stn_min, north_min, imask)
    call imap%stack_on_peaks(sts_max, south_max, imask)
    call imap%stack_on_peaks(sts_min, south_min, imask)

    call north_max%get_all_radial_profiles()
    call north_min%get_all_radial_profiles()
    call south_max%get_all_radial_profiles()
    call south_min%get_all_radial_profiles()
    
    north_max%fr = north_max%fr*1.d6
    north_min%fr = north_min%fr*1.d6
    south_max%fr = south_max%fr*1.d6
    south_min%fr = south_min%fr*1.d6    
  end subroutine stack_imap

  
  subroutine load_imap(i)
    COOP_INT i
    logical,save::first_load = .true.
    if(iloaded) return
    if(i.eq.0)then
       call imap%read(imap_file, nmaps_wanted = sto_max%nmaps, nmaps_to_read = 1)
    else
       if(first_load)then
          call imap%read(trim(mapdir)//"dx11_v2_"//trim(cc_method)//"_int_cmb_mc_"//trim(coop_Ndigits(i-1, 5))//trim(postfix), nmaps_to_read = 1, nmaps_wanted = sto_max%nmaps)
          call inoise%read(trim(mapdir)//"dx11_v2_"//trim(cc_method)//"_int_noise_mc_"//trim(coop_Ndigits(i-1, 5))//trim(postfix), nmaps_to_read = 1, nmaps_wanted = sto_max%nmaps)
       else          
          call imap%read(trim(mapdir)//"dx11_v2_"//trim(cc_method)//"_int_cmb_mc_"//trim(coop_Ndigits(i-1, 5))//trim(postfix), nmaps_to_read = 1, known_size = .true.)
          call inoise%read(trim(mapdir)//"dx11_v2_"//trim(cc_method)//"_int_noise_mc_"//trim(coop_Ndigits(i-1, 5))//trim(postfix), nmaps_to_read = 1, known_size = .true.)
          first_load = .false.
       endif
       imap%map(:, 1) = imap%map(:, 1) + inoise%map(:, 1)       
    endif
    if(imap%nmaps .eq. 3)then
       imap%map(:,1) = imap%map(:,1)*imask%map(:,1)
       call imap%iqu2TQTUT()
    endif
  end subroutine load_imap
  
#else
  print*, "You need to install healpix"
#endif  
end program Exp_spots
