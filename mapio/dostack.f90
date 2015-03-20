program Stacking_Maps
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX
  logical::use_mask = .true.
  logical::remove_mono = .false.
  
  COOP_STRING::stack_field_name = "T"
  COOP_STRING::map_file = "act15/act15_i_hp_230_270_smoothed_fwhm5arcmin.fits"
  COOP_STRING::imask_file = "act15/act15_imask.fits"
  COOP_STRING::polmask_file = "act15/act15_polmask.fits"
  COOP_STRING::peak_file = "peaks/act15_nu1_fp.dat"
  COOP_UNKNOWN_STRING,parameter::mask_file_force_to_use = ""
  
  COOP_INT,parameter::n = 36
  COOP_REAL,parameter::r_degree  = 2.d0
  COOP_REAL,parameter::dr = 2.d0*sin(r_degree*coop_SI_degree/2.d0)/n
  logical::makepdf = .false.
  type(coop_stacking_options)::sto
  type(coop_healpix_patch)::patch
  type(coop_healpix_maps)::hgm, mask, pmap
  COOP_STRING::output 
  COOP_INT i, m
  COOP_REAL::zmin1 = 1.1e31
  COOP_REAL::zmax1 = -1.1e31
  COOP_REAL::zmin2 = 1.1e31
  COOP_REAL::zmax2 = -1.1e31
  COOP_STRING::line  
  type(coop_asy)::fig
  COOP_REAL::tmax
  
  output = "stacked/act15_"//trim(stack_field_name)//"_nu1_5arcmin_fpts"
  if(iargc() .ge. 6)then
     use_mask = .true.
     map_file = coop_InputArgs(1)
     imask_file = coop_InputArgs(2)
     polmask_file = coop_InputArgs(3)     
     peak_file = coop_InputArgs(4)
     stack_field_name = coop_InputArgs(5)
     output = coop_InputArgs(6)
     if(iargc() .ge. 8)then
        line = coop_InputArgs(7)
        read(line, *) zmin1
        line = coop_InputArgs(8)
        read(line, *) zmax1
     endif
     if(iargc() .ge. 10)then
        line = coop_InputArgs(9)
        read(line, *) zmin2
        line = coop_InputArgs(10)
        read(line, *) zmax2
     endif
  else
     makepdf = .true.
     coop_healpix_patch_default_want_caption = .true.
     coop_healpix_patch_default_want_label  = .true.
     coop_healpix_patch_default_figure_width = 3.5
     coop_healpix_patch_default_figure_height = 3.   
     coop_healpix_mask_tol = 0.d0
     if(.not. use_mask)then
        write(*,*) "Warning: not using the mask"
     endif
  endif
  call sto%import(peak_file)
  call hgm%read(map_file)
  call patch%init(stack_field_name, n, dr)
  if(use_mask)then
     if(patch%tbs%mask_int .and. .not. patch%tbs%mask_pol)then
        call mask%read(imask_file, nmaps_wanted = 1, spin = (/ 0 /))
     elseif(patch%tbs%mask_pol .and. .not. patch%tbs%mask_int)then
        call mask%read(polmask_file, nmaps_wanted = 1, spin = (/ 0 /))
     elseif(trim(mask_file_force_to_use) .ne. "")then
        call mask%read(mask_file_force_to_use, nmaps_wanted = 1, spin = (/ 0 /))
     else
        stop "For unclassified stacking maps you have to specify the mask explicitly"
     endif
     if(mask%nside .ne. hgm%nside) stop "mask and map must have the same nside"
     do i=1, hgm%nmaps
        hgm%map(:, i) = hgm%map(:, i)*mask%map(:, 1)
        if(remove_mono) hgm%map(:, i) = hgm%map(:, i) - sum(dble(hgm%map(:,i)))/sum(dble(mask%map(:,1))) !!remove monopole
     enddo
  endif
  tmax = maxval(hgm%map(:,1))
  if(tmax .lt. 1.d0 .and. tmax .gt. 1.d-5)then
     hgm%map = hgm%map*1.e6
  endif
  print*, "stacking on "//COOP_STR_OF(sto%peak_pix%n)//" peaks"  
  if(use_mask)then
     call hgm%stack_on_peaks(sto, patch, mask)
  else
     call hgm%stack_on_peaks(sto, patch)
  endif
  patch%caption = "stacked on "//COOP_STR_OF(sto%peak_pix%n)//" "//trim(sto%caption)
  select case(patch%nmaps)
  case(1)
     patch%tbs%zmin(1) = zmin1
     patch%tbs%zmax(1) = zmax1     
     call patch%plot(1, trim(adjustl(output))//".txt")
     if(makepdf)call system("../utils/fasy.sh "//trim(adjustl(output))//".txt")
  case default
     patch%tbs%zmin(1) = zmin1
     patch%tbs%zmax(1) = zmax1
     patch%tbs%zmin(2) = zmin2          
     patch%tbs%zmax(2) = zmax2          
     do i=1, patch%nmaps
        call patch%plot(i, trim(adjustl(output))//"_"//COOP_STR_OF(i)//".txt")
        if(makepdf)call system("../utils/fasy.sh "//trim(adjustl(output))//"_"//COOP_STR_OF(i)//".txt")        
     enddo
  end select
  call patch%get_all_radial_profiles()
  select case(patch%nmaps)
  case(1)
     do m = 0, patch%mmax, 2
        call fig%open(trim(adjustl(output))//"_m"//COOP_STR_OF(m)//".txt")
        call fig%init(xlabel="$r$", ylabel="radial profile")
        call coop_asy_curve(fig, patch%r, patch%fr(:, m/2, 1))
        call fig%close()
        call fig%open(trim(adjustl(output))//"_m"//COOP_STR_OF(m)//".dat")
        do i=0, patch%n
           write(fig%unit, "(2E14.5)") patch%r(i), patch%fr(i, m/2, 1)
        enddo
        call fig%close()        
     enddo
  case(2)
     do m = 0, patch%mmax, 2
        call fig%open(trim(adjustl(output))//"_m"//COOP_STR_OF(m)//".txt")
        call fig%init(xlabel="$r$", ylabel="radial profile")
        if(m.ne.0)then
           call coop_asy_curve(fig, patch%r, (patch%fr(:, m/2, 1)+patch%fr(:, m/2, 2))/2.d0)
        else
           call coop_asy_curve(fig, patch%r, patch%fr(:, m/2, 1) )
        endif
        call fig%close()
        call fig%open(trim(adjustl(output))//"_m"//COOP_STR_OF(m)//".dat")
        do i=0, patch%n
           write(fig%unit, "(3E14.5)") patch%r(i), patch%fr(i, m/2, :)
        enddo
        call fig%close()        
     enddo
  end select
#else
  print*, "You need to install healpix"
#endif  
end program Stacking_Maps
