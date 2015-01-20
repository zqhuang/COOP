program Stacking_Maps
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX
  logical::do_max = .true.
  COOP_STRING::stack_field_name = "QU"
  COOP_STRING::map_file = "planck14/dx11_v2_smica_pol_case1_cmb_010a_1024.fits"
  COOP_STRING::peak_file = "peaks/smica_fwhm10_Tmax_ORIENTNULL_nu0.dat"
  COOP_STRING::imask_file = "planck14/dx11_v2_common_int_mask_010a_1024.fits"
  COOP_STRING::polmask_file = "planck14/dx11_v2_common_pol_mask_010a_1024.fits"
  COOP_INT,parameter::n = 36
  COOP_REAL,parameter::r_degree  = 2.d0
  COOP_REAL,parameter::dr = 2.d0*sin(r_degree*coop_SI_degree/2.d0)/n
  
  type(coop_stacking_options)::sto
  type(coop_healpix_patch)::patch
  type(coop_healpix_maps)::hgm, mask, pmap
  COOP_STRING::output = "stacked/sample"
  COOP_INT i, m
  COOP_REAL::zmin1 = 1.1e31
  COOP_REAL::zmax1 = -1.1e31
  COOP_REAL::zmin2 = 1.1e31
  COOP_REAL::zmax2 = -1.1e31
  COOP_STRING :: line
  
  type(coop_asy)::fig
  if(iargc() .ge. 6)then
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
        
  endif
  call sto%import(peak_file)
  call hgm%read(map_file)
  call patch%init(stack_field_name, n, dr)
  if(all(patch%tbs%spin.eq.0))then
     call mask%read(imask_file)
  else
     call mask%read(polmask_file)
  endif
  if(mask%nside .ne. hgm%nside) stop "mask and map must have the same nside"
  do i=1, hgm%nmaps
     hgm%map(:, i) = hgm%map(:,i)*mask%map(:, 1)
  enddo
  if(maxval(hgm%map(:,1)) .lt. 1.d-2)then
     hgm%map = hgm%map*1.e6
  endif
  print*, "stacking on "//COOP_STR_OF(sto%peak_pix%n)//" peaks"  
  
  call hgm%stack_on_peaks(sto, patch, mask)
  patch%caption = "stacked on "//COOP_STR_OF(sto%peak_pix%n)//" "//trim(sto%caption)
  select case(patch%nmaps)
  case(1)
     patch%tbs%zmin(1) = zmin1
     patch%tbs%zmax(1) = zmax1     
     call patch%plot(1, trim(adjustl(output))//".txt")
  case default
     patch%tbs%zmin(1) = zmin1
     patch%tbs%zmax(1) = zmax1
     patch%tbs%zmin(2) = zmin2          
     patch%tbs%zmax(2) = zmax2          
     do i=1, patch%nmaps
        call patch%plot(i, trim(adjustl(output))//"_"//COOP_STR_OF(i)//".txt")
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
        call coop_asy_curve(fig, patch%r, patch%fr(:, m/2, 1)+patch%fr(:, m/2, 2))
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
