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


  COOP_REAL, parameter::pre_smooth = 30*coop_SI_arcmin
  COOP_UNKNOWN_STRING, parameter::color_table = "Rainbow"
  COOP_UNKNOWN_STRING, parameter::spot_type = "Tmax"
  COOP_UNKNOWN_STRING, parameter::stack_type = "T"
  COOP_REAL, parameter::threshold = 0
  COOP_INT, parameter::mmax = 0
  COOP_INT, parameter::n = 20
  COOP_REAL, parameter::dr = 15.*coop_SI_arcmin
  COOP_UNKNOWN_STRING, parameter::map_file = "commander/commander_dx11d2_temp_n2048_fullres_hybrid_v3_full_cmb.fits"
  COOP_UNKNOWN_STRING, parameter::spots_mask_file  = "commander/commander_dx11d2_mask_temp_n2048_fullres_v3.fits"
  COOP_UNKNOWN_STRING, parameter::stack_mask_file  = "commander/commander_dx11d2_mask_temp_n2048_fullres_v3.fits"
  COOP_UNKNOWN_STRING, parameter::prefix = "hsloutput/"
  type(coop_healpix_patch)::patch_s, patch_n
  integer,parameter::scan_nside = 4, n_sim = 3, imap = 1,  nmaps_wanted = 1
  integer run_id, i, nlines, n_larger_chisq, weight
  type(coop_healpix_maps)::map, spots_mask, stack_mask
  COOP_REAL diff(0:n), chisq, prob, hdir(2), kdata, bdata, kmean, bmean, cov(2,2)
  COOP_STRING::fmt, fr_file, log_file, fig_file
  type(coop_list_integer)::listpix
  type(coop_list_real)::listangle
  type(coop_file) fp
  type(coop_asy) fig
  COOP_REAL,dimension(n_sim)::ksim, bsim

  run_id = coop_str2int(coop_InputArgs(1))
  if(run_id .ge.  scan_nside**2*12)then
     write(*,*) "run id must not exceed ", scan_nside**2*12 - 1
     stop
  endif
  call pix2ang_ring(scan_nside, run_id, hdir(1), hdir(2))

  !!read mask and map
  call spots_mask%read(spots_mask_file, nmaps_wanted = 1)
  call stack_mask%read(stack_mask_file, nmaps_wanted = 1)

  call patch_n%init(stack_type, n, dr, mmax = mmax)
  patch_s = patch_n
  fr_file = prefix//"fr_"//COOP_STR_OF(scan_nside)//"_"//COOP_STR_OF(run_id)//".txt"
  fig_file = prefix//"fig_"//COOP_STR_OF(scan_nside)//"_"//COOP_STR_OF(run_id)//".txt"
  log_file = prefix//"log_"//COOP_STR_OF(scan_nside)//"_"//COOP_STR_OF(run_id)//".txt"
  fmt = "("//trim(coop_num2str(n+1))//"G16.7)"
  weight = 0
  if(coop_file_exists(trim(fr_file)))then
     nlines = coop_file_numlines(trim(fr_file))
     if(nlines .gt. 0)then
        call fp%open(trim(fr_file), "r")
        do i = 1, nlines
           read(fp%unit, trim(fmt)) diff
           weight = weight + 1
           call coop_linear_least_square_fit(n+1, patch_s%r, diff, ksim(weight), bsim(weight))
           if(weight .ge. n_sim) exit
        enddo
        call fp%close()
        write(*,*) "Loaded "//trim(coop_num2str(weight))//" lines from checkpoint"
     endif
  endif


  call fp%open(trim(fr_file), "a")
  do while(weight .lt. n_sim)
     weight = weight + 1
     call map%read(trim(sim_file_name(weight)), nmaps_wanted = nmaps_wanted)
     call map%mask(spots_mask)
     call map%smooth(pre_smooth)
     call map%get_listpix(listpix, listangle, spot_type, threshold, spots_mask)
     call map%convert2ring()
     call map%stack_north_south(patch_n, patch_s, listpix, listangle, hdir, stack_mask)
     call patch_n%get_radial_profile(imap = imap, m = 0)
     call patch_s%get_radial_profile(imap = imap, m = 0)
     diff = patch_n%fr(:, 0, imap) - patch_s%fr(:, 0, imap)
     call coop_linear_least_square_fit(n+1, patch_s%r, diff, ksim(weight), bsim(weight))
     write(fp%unit, trim(fmt)) diff
     flush(fp%unit)
  enddo
  call fp%close()

  call map%read(map_file)
  call map%mask(spots_mask)
  call map%smooth(pre_smooth)
  call map%get_listpix(listpix, listangle, spot_type, threshold, spots_mask)
  call map%convert2ring()
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
  call fig%init(ylabel = "$\delta T (\mu K) $", xlabel = "$dT /dr (\mu K/{\rm rad})$")
  call coop_asy_dots(fig, ksim, bsim, "black")
  call coop_asy_dot(fig, kdata, bdata, "red", "x")
  call coop_asy_label(fig, "x : data",  fig%xmin+(fig%xmax - fig%xmin)*0.12, fig%ymin + 0.86*(fig%ymax - fig%ymin), color = "red")
  call coop_asy_label(fig, "$\bullet$ : sims",  fig%xmin+(fig%xmax - fig%xmin)*0.15, fig%ymin + 0.82*(fig%ymax - fig%ymin), color = "black")
  call coop_asy_label(fig, "direction: $l ="//trim(coop_num2str(nint(hdir(2)/coop_SI_degree)))//"^{\circ},  b="//trim(coop_num2str(nint((coop_pio2-hdir(1))/coop_SI_degree)))//"^{\circ}$",  fig%xmin+(fig%xmax - fig%xmin)*0.2, fig%ymin + 0.78*(fig%ymax - fig%ymin), color = "black")
  call fig%close()

contains

  function sim_file_name(i)
    COOP_INT i
    COOP_STRING sim_file_name
    select case(i)
    case(1)
       sim_file_name = "ffp7/ffp7_smica_cmb_0001_2048_debeam.fits"
    case(2)
       sim_file_name = "planck/smica_inp_cmb.fits"
    case(3)
       sim_file_name = "simu/simulate_iqu_n2048.fits"
    end select
  end function sim_file_name

end program test
