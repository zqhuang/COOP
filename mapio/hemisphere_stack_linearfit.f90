program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  implicit none
#include "constants.h"
  COOP_REAL, parameter::pre_smooth = 15*coop_SI_arcmin
  COOP_UNKNOWN_STRING, parameter::prefix = "predx11"
  COOP_UNKNOWN_STRING, parameter::color_table = "Planck"
  COOP_UNKNOWN_STRING, parameter::spot_type = "Tmax"
  COOP_UNKNOWN_STRING, parameter::stack_type = "T"
  COOP_REAL, parameter::threshold = 0
  COOP_INT, parameter::mmax = 0
  COOP_INT, parameter::n = 30
  COOP_REAL, parameter::dr = 10.*coop_SI_arcmin
  COOP_INT, parameter:: imap = 1
  integer,parameter::n_sim = 1000
  COOP_STRING::fmt, fmtscreen
  COOP_UNKNOWN_STRING, parameter::resol = "1024"
  COOP_UNKNOWN_STRING, parameter::map_file = prefix//"/"//prefix//"_iqu"//resol//".fits"
  COOP_UNKNOWN_STRING, parameter::imask_file  = prefix//"/"//prefix//"_imask"//resol//".fits"
  COOP_UNKNOWN_STRING, parameter::polmask_file  = prefix//"/"//prefix//"_polmask"//resol//".fits"
  COOP_UNKNOWN_STRING, parameter::fr_file = prefix//"_"//stack_type//"_sim_Fr"//resol//".dat"
  type(coop_healpix_patch)::patch_s, patch_n
  type(coop_healpix_maps)::map, sim, tmp
  type(coop_healpix_maps)::imask, polmask
  type(coop_healpix_maps)::stack_mask, spots_mask
  COOP_REAL diff(0:n), chisq,  hdir(2), kdata, bdata, kmean, bmean, ksig, bsig
  COOP_REAL,dimension(n_sim)::kdiff, bdiff
  COOP_INT :: nspots, nlines, il, i, j, step, weight
  type(coop_asy)::fig
  type(coop_file)::fp
  type(coop_list_integer)::listpix
  type(coop_list_real)::listangle
  type(coop_file) fpsim
  integer,parameter::ncols_used = 30
  COOP_REAL, parameter::epsilon = 1.d-6
  COOP_INT nmaps_wanted, ismall, ilarge

  !!read mask and map
  call imask%read(imask_file, nmaps_wanted = 1)
  call polmask%read(polmask_file, nmaps_wanted = 1)
  nmaps_wanted = 1
  select case(trim(stack_type))
  case("I","T")
     stack_mask = imask
  case default
     stack_mask = polmask
     nmaps_wanted = 3
  end select
  select case(trim(spot_type))
  case("Tmax", "Tmin", "Tmax_QTUTOrient", "Tmin_QTUTOrient")
     spots_mask = imask
  case default
     spots_mask = polmask
     nmaps_wanted = 3
  end select
  call map%read(map_file, nmaps_wanted = nmaps_wanted)
  if(nmaps_wanted .eq. 3)then
     call map%mask(imask, polmask)
  else
     call map%mask(imask)
  endif
  if(pre_smooth .gt. coop_SI_arcmin)then
     call map%smooth(pre_smooth)
  endif
  if(nmaps_wanted .eq. 3)then
     call map%get_fullCls(imask, polmask)
  else
     call map%get_fullCls(imask)
  endif
  call coop_healpix_lb2ang(l_deg = 226.d0, b_deg = -17.d0, theta = hdir(1), phi = hdir(2))

  sim = map
  tmp = sim

  call coop_random_init()
  call patch_n%init(stack_type, n, dr, mmax = mmax)
  patch_s = patch_n

  fmt = "("//trim(coop_num2str(n+1))//"G16.7)"
  fmtscreen = "("//trim(coop_num2str(n+1))//"F7.2)"
  weight = 0
  if(coop_file_exists(fr_file))then
     nlines = coop_file_numlines(fr_file) 
     if(nlines .gt. 0)then
        call fpsim%open(fr_file, "r")
        do il=1, nlines
           read(fpsim%unit, trim(fmt)) diff
           weight = weight + 1
           call coop_linear_least_square_fit(n+1, patch_s%r, diff, kdiff(weight), bdiff(weight))
           if(weight .ge. n_sim) exit
        enddo
        call fpsim%close()
        write(*,*) "Loaded "//trim(coop_num2str(weight))//" lines from checkpoint"
     endif
  endif


  call fpsim%open(fr_file, "a")
  do while(weight .lt. n_sim)
     if(mod(weight, 20) .eq. 0)then
        if(nmaps_wanted .eq. 3)then
           call map%get_fullCls(imask, polmask)
        else
           call map%get_fullCls(imask)
        endif
        sim%cl = map%cl
     endif
     call sim%simulate()
     tmp%map = sim%map
     tmp%ordering = sim%ordering
     weight = weight + 1
     call tmp%get_listpix(listpix, listangle, spot_type, threshold, spots_mask)
     call sim%stack_north_south(patch_n, patch_s, listpix, listangle, hdir, stack_mask)
     write(*,"(A)") "sim #"//trim(coop_num2str(weight))//", num_N = "//trim(coop_num2str(patch_n%nstack_raw))//", num_S = "//trim(coop_num2str(patch_s%nstack_raw))
     call patch_n%get_radial_profile(imap = imap, m = 0)
     call patch_s%get_radial_profile(imap = imap, m = 0)
     diff = patch_n%fr(:, 0, imap) - patch_s%fr(:, 0, imap)
     call coop_linear_least_square_fit(n+1, patch_s%r, diff, kdiff(weight), bdiff(weight))
     write(fpsim%unit, trim(fmt)) diff
     flush(fpsim%unit)
     write(*, trim(fmtscreen)) diff
  enddo
  call fpsim%close()

  sim = map
  call sim%get_listpix(listpix, listangle, spot_type, threshold, spots_mask)
  call map%stack_north_south(patch_n, patch_s, listpix, listangle, hdir, stack_mask )
  call patch_s%get_radial_profile(imap = imap, m = 0)
  call patch_n%get_radial_profile(imap = imap, m = 0)
  diff = patch_n%fr(:, 0, imap) - patch_s%fr(:, 0, imap)
  call coop_linear_least_square_fit(n+1, patch_s%r, diff, kdata, bdata)
  kmean = sum(kdiff)/n_sim
  bmean  = sum(bdiff)/n_sim
  ksig = sqrt( sum((kdiff - kmean)**2)/n_sim )
  bsig = sqrt( sum((bdiff - bmean)**2)/n_sim )
  print*, (kdata - kmean)/ksig, count(kdiff .gt. kdata)/real(n_sim)
  print*, (bdata-bmean)/bsig, count(bdiff .gt. bdata)/real(n_sim)

end program test
