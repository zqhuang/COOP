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
  COOP_INT, parameter::dir_l = 226
  COOP_INT, parameter::dir_b = -17
  COOP_UNKNOWN_STRING, parameter::prefix="simu"
  COOP_UNKNOWN_STRING, parameter::color_table = "Rainbow"
  COOP_UNKNOWN_STRING, parameter::spot_type = "Tmax"
  COOP_UNKNOWN_STRING, parameter::stack_type = "T"
  COOP_REAL, parameter::threshold = 0
  COOP_INT, parameter::mmax = 0
  COOP_INT, parameter::n = 30
  COOP_REAL, parameter::dr = 10.*coop_SI_arcmin
  COOP_INT, parameter:: imap = 1
  integer,parameter::n_sim = 420
  COOP_STRING::fmt, fmtscreen
  COOP_UNKNOWN_STRING, parameter::resol = "256"
  COOP_UNKNOWN_STRING, parameter::map_file = prefix//"/"//prefix//"_iqu"//resol//".fits"
  COOP_UNKNOWN_STRING, parameter::imask_file  = prefix//"/"//prefix//"_imask"//resol//".fits"
  COOP_UNKNOWN_STRING, parameter::polmask_file  = prefix//"/"//prefix//"_polmask"//resol//".fits"
  COOP_STRING::fr_file, log_file, fig_file
  type(coop_healpix_patch)::patch_s, patch_n
  type(coop_healpix_maps)::map, sim, tmp
  type(coop_healpix_maps)::imask, polmask
  type(coop_healpix_maps)::stack_mask, spots_mask
  COOP_REAL diff(0:n), chisq, prob, hdir(2), kdata, bdata, kmean, bmean, cov(2,2)
  COOP_REAL,dimension(n_sim)::ksim, bsim
  COOP_INT :: nspots, nlines, il, i, j, weight, n_larger_chisq
  type(coop_asy)::fig
  type(coop_file)::fp
  type(coop_list_integer)::listpix
  type(coop_list_real)::listangle
  type(coop_file) fpsim
  integer::scan_nside
  COOP_INT nmaps_wanted
  COOP_SHORT_STRING:: input, input2
  integer run_id
  input = trim(coop_inputArgs(1))
  if(trim(input) .eq. "") then
     scan_nside = 0
     run_id = 0
     call coop_healpix_lb2ang(l_deg = dble(dir_l), b_deg = dble(dir_b), theta = hdir(1), phi = hdir(2))
     fr_file = prefix//"/"//prefix//resol//"_"//stack_type//"_l"//trim(coop_num2str(dir_l))//"_b"//trim(coop_num2str(dir_b))//".txt"
     log_file = prefix//"/"//prefix//resol//"_"//stack_type//"_l"//trim(coop_num2str(dir_l))//"_b"//trim(coop_num2str(dir_b))//".log"
     fig_file = prefix//"/figure_"//prefix//resol//"_"//stack_type//"_l"//trim(coop_num2str(dir_l))//"_b"//trim(coop_num2str(dir_b))//".txt"
  else
     input2=trim(coop_inputArgs(2))
     if(trim(input2).eq."")then
        scan_nside = 2
     else
        read(input2, *) scan_nside
     endif
     read(input,*) run_id
     write(*,*) "nside, id = ", scan_nside, run_id
     if(run_id .ge.  scan_nside**2*12)then
        write(*,*) "run id must not exceed ", scan_nside**2*12 - 1
        stop
     endif
     call pix2ang_ring(scan_nside, run_id, hdir(1), hdir(2))
     fr_file = prefix//"/"//prefix//resol//"_"//stack_type//"_"//trim(input2)//"id"//trim(input)//".txt"
     log_file = prefix//"/"//prefix//resol//"_"//stack_type//"_"//trim(input2)//"id"//trim(input)//".log"
     fig_file = prefix//"/figure_"//prefix//resol//"_"//stack_type//"_"//trim(input2)//"id"//trim(input)//".txt"
  endif
  write(*, "(A,2E14.3)") "theta, phi = ", hdir

  !!read mask and map
  call imask%read(imask_file, nmaps_wanted = 1)
  nmaps_wanted = 1
  select case(trim(stack_type))
  case("I","T")
     stack_mask = imask
  case default
     nmaps_wanted = 3
     call polmask%read(polmask_file, nmaps_wanted = 1)
     stack_mask = polmask
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
  


  sim = map
  tmp = sim

  call coop_random_init()
  call patch_n%init(stack_type, n, dr, mmax = mmax)
  patch_s = patch_n

  fmt = "("//trim(coop_num2str(n+1))//"G16.7)"
  fmtscreen = "("//trim(coop_num2str(n+1))//"F7.2)"
  weight = 0
  if(coop_file_exists(trim(fr_file)))then
     nlines = coop_file_numlines(trim(fr_file))
     if(nlines .gt. 0)then
        call fpsim%open(trim(fr_file), "r")
        do il=1, nlines
           read(fpsim%unit, trim(fmt)) diff
           weight = weight + 1
           call coop_linear_least_square_fit(n+1, patch_s%r, diff, ksim(weight), bsim(weight))
           if(weight .ge. n_sim) exit
        enddo
        call fpsim%close()
        write(*,*) "Loaded "//trim(coop_num2str(weight))//" lines from checkpoint"
     endif
  endif


  call fpsim%open(trim(fr_file), "a")
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
     call coop_linear_least_square_fit(n+1, patch_s%r, diff, ksim(weight), bsim(weight))
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
  kmean = sum(ksim)/n_sim
  bmean  = sum(bsim)/n_sim
  cov(1,1) = sum((ksim-kmean)**2)/n_sim
  cov(2,2) = sum((bsim-bmean)**2)/n_sim
  cov(1,2) = sum((ksim-kmean)*(bsim-bmean))/n_sim
  cov(2,1) = cov(1,2)
  call coop_matsym_inverse_small(2, cov)
  chisq = (kdata - kmean)**2*cov(1,1) + (bdata - bmean)**2*cov(2,2) + (kdata -kmean)*(bdata-bmean)*2.d0*cov(1,2)
  n_larger_chisq = 0
  do i=1, n_sim
     if((ksim(i) - kmean)**2*cov(1,1) + (bsim(i) - bmean)**2*cov(2,2) + (ksim(i) -kmean)*(bsim(i)-bmean)*2.d0*cov(1,2) .gt. chisq) then
        n_larger_chisq = n_larger_chisq + 1
     endif
  enddo
  write(*,"(A)") "chi^2 = "//trim(coop_num2str(chisq, "(F11.3)"))
  prob = dble(n_larger_chisq)/n_sim
  write(*,"(A)") "The probability of seeing a larger chi^2 in simulations is: "//trim(coop_num2str(100.*prob,"(F10.2)"))//"%"
  

  call fp%open(trim(log_file), "w")
  write(fp%unit, "(2I6, 4E16.7)") scan_nside, run_id, hdir, chisq, prob
  call fp%close()
     
  call fig%open(trim(fig_file))
  call fig%init(ylabel = "$\delta T (\mu K) $", xlabel = "$dT /dr (\mu K/{\rm rad})$", caption = "black dots are simulations, red x is data")
  call coop_asy_dots(fig, ksim, bsim, "black")
  call coop_asy_dot(fig, kdata, bdata, "red", "x")
!!  call coop_asy_label(fig, "Map: "//prefix,  fig%xmin+(fig%xmax - fig%xmin)*0.12, fig%ymin + 0.9*(fig%ymax - fig%ymin))
  call coop_asy_label(fig, "x : mock data",  fig%xmin+(fig%xmax - fig%xmin)*0.12, fig%ymin + 0.86*(fig%ymax - fig%ymin), color = "red")
  call coop_asy_label(fig, "$\bullet$ :  isotropized sims",  fig%xmin+(fig%xmax - fig%xmin)*0.15, fig%ymin + 0.82*(fig%ymax - fig%ymin), color = "black")
  call coop_asy_label(fig, "direction: $l ="//trim(coop_num2str(nint(hdir(2)/coop_SI_degree)))//"^{\circ},  b="//trim(coop_num2str(nint((coop_pio2-hdir(1))/coop_SI_degree)))//"^{\circ}$",  fig%xmin+(fig%xmax - fig%xmin)*0.2, fig%ymin + 0.78*(fig%ymax - fig%ymin), color = "black")
  call fig%close()

end program test
