program Exp_spots
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX

  COOP_STRING::cc_method 
  
  COOP_STRING::output, line
  COOP_REAL::threshold

  COOP_INT::n_sim = 30
  COOP_STRING,parameter::peak_name = "$T$"
  COOP_STRING::orient_name 
  COOP_STRING::stack_field_name = "T"
  COOP_INT,parameter::n = 36
  COOP_REAL,parameter::r_degree  = 2.d0
  COOP_REAL,parameter::dr = 2.d0*sin(r_degree*coop_SI_degree/2.d0)/n

  COOP_REAL::r(0:n), pfr(0:n), pfr0(0:n)
  COOP_INT::count_p
  COOP_UNKNOWN_STRING,parameter::outputdir = "st/"  
  COOP_UNKNOWN_STRING,parameter::mapdir = "massffp8/"
  COOP_UNKNOWN_STRING,parameter::postfix = "_020a_0512.fits"
  
  COOP_STRING::imap_file, polmap_file, imask_file, polmask_file
  
  type(coop_stacking_options)::sto_max, sto_min
  type(coop_healpix_patch)::patch_max, patch_min
  type(coop_healpix_maps)::imap, imask, polmask, inoise, polnoise, polmap
  logical::iloaded = .false.
  logical::polloaded  = .false.
  type(coop_file)::fp
  type(coop_asy)::fig
  COOP_INT i
  COOP_INT ind, ind_done
  
  if(iargc() .lt. 3)then
     write(*,*) "Syntax:"
     write(*,*) "./SST cc_method nu n_sim [Orient] [cut_cold_spot]"
     write(*,*) "Examples:"     
     write(*,*) "./SST smica 0.5 1000"
     write(*,*) "./SST nilc 0. 100 T T"
     write(*,*) "./SST commander 2. 100 F T"          
     stop
  endif
  cc_method = trim(coop_inputArgs(1))
  imap_file = "planck14/dx11_v2_"//trim(cc_method)//"_int_cmb"//postfix
  polmap_file = "planck14/dx11_v2_"//trim(cc_method)//"_pol_case1_cmb_hp_20_40"//postfix
  
  line = coop_inputArgs(2)
  read(line, *) threshold
  line = coop_inputArgs(3)
  read(line,*) n_sim

  if(trim(coop_inputArgs(4)).eq."T")then
     orient_name = "$(Q_T, U_T)$"
  else
     orient_name = "NULL"
  endif
  
  if(trim(coop_inputArgs(5)).eq."T")then
     imask_file = "planck14/coldspot_mask"//postfix
     output = trim(outputdir)//trim(cc_method)//"_nu"//trim(coop_num2goodstr(threshold,"-","pt"))//"_"//trim(coop_str_numalpha(peak_name))//"_Orient"//trim(coop_str_numalpha(orient_name))//"_CutColdSpot.dat"            
  else
     imask_file = "planck14/dx11_v2_common_int_mask"//postfix
     output = trim(outputdir)//trim(cc_method)//"_nu"//trim(coop_num2goodstr(threshold,"-","pt"))//"_"//trim(coop_str_numalpha(peak_name))//"_Orient"//trim(coop_str_numalpha(orient_name))//".dat"                 
  endif
  polmask_file = "planck14/dx11_v2_common_pol_mask"//postfix
  
  


  
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
    
  call patch_max%init(trim(stack_field_name), n, dr)
  patch_min = patch_max
  
  ind_done = -1  
  call imask%read(imask_file, nmaps_wanted = 1, spin = (/ 0 /) )

  if(coop_file_exists(output))then
     call fp%open(output, "ru")
     do i = 0, n_sim
        read(fp%unit, ERR = 100, END = 100) ind, patch_max%fr, patch_min%fr
        ind_done = ind
     enddo
100  call fp%close()
  endif
  
  call fp%open(output, "u")
  call fig%open("fr_diff.txt")
  call fig%init(xlabel = "$\varpi$", ylabel = "$T_0^{\rm hot} + T_0^{\rm cold}$")
  if(ind_done .ge. 0)then
     print*, "loaded "//COOP_STR_OF(ind_done+1)//" stacked maps"
  endif
  count_p = 0
  do ind = 0, n_sim
     if(ind.gt.ind_done)then
        print*, "stacking map#"//COOP_STR_OF(ind)                
        call load_imap(ind)
        call stack_imap()

        write(fp%unit) ind, patch_max%fr, patch_min%fr
     else
        read(fp%unit) i, patch_max%fr, patch_min%fr
     endif
     pfr = (patch_min%fr(:, 0, 1)) !-patch_min%fr(:,0,1))/2.d0
     if(ind.eq.0)then
        pfr0 = pfr
!        call fig%curve(r, pfr, color = "red", linetype= "solid", linewidth = 1.5)
     else
        call fig%curve(r, pfr - pfr0, color = "blue", linetype= "dotted", linewidth = 1.)
        if(sum(pfr-pfr0).lt.0.d0)then
           count_p = count_p + 1
        endif
     endif
  enddo
  print*, COOP_STR_OF(count_p)//" in "//COOP_STR_OF(n_sim)//" has negative S"
  call fig%close()
  call fp%close()
  

contains

  subroutine stack_imap()
    call imap%get_peaks(sto_max, mask = imask, restore = .false.)
    call imap%get_peaks(sto_min, mask = imask, restore = .false.)
    call imap%stack_on_peaks(sto_max, patch_max, imask)
    call imap%stack_on_peaks(sto_min, patch_min, imask)
    call patch_max%get_all_radial_profiles()
    call patch_min%get_all_radial_profiles()
    patch_max%fr = patch_max%fr*1.d6
    patch_min%fr = patch_min%fr*1.d6    
  end subroutine stack_imap

  subroutine load_imap(i)
    COOP_INT i
    logical,save::first_load = .true.
    if(iloaded) return
    if(i.eq.0)then
       call imap%read(imap_file, nmaps_wanted = sto_max%nmaps, nmaps_to_read = 1)
    else
       if(first_load)then
          call imap%read(trim(mapdir)//"cmb/dx11_v2_"//trim(cc_method)//"_int_cmb_mc_"//trim(coop_Ndigits(i-1, 5))//trim(postfix), nmaps_to_read = 1, nmaps_wanted = sto_max%nmaps)
          call inoise%read(trim(mapdir)//"noise/dx11_v2_"//trim(cc_method)//"_int_noise_mc_"//trim(coop_Ndigits(i-1, 5))//trim(postfix), nmaps_to_read = 1, nmaps_wanted = sto_max%nmaps)
       else          
          call imap%read(trim(mapdir)//"cmb/dx11_v2_"//trim(cc_method)//"_int_cmb_mc_"//trim(coop_Ndigits(i-1, 5))//trim(postfix), nmaps_to_read = 1, known_size = .true.)
          call inoise%read(trim(mapdir)//"noise/dx11_v2_"//trim(cc_method)//"_int_noise_mc_"//trim(coop_Ndigits(i-1, 5))//trim(postfix), nmaps_to_read = 1, known_size = .true.)
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
