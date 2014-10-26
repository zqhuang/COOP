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
  COOP_INT, parameter::n_sim = 1000
  COOP_UNKNOWN_STRING, parameter::color_table = "Rainbow"
  COOP_UNKNOWN_STRING, parameter::spot_type = "Tmax"
  COOP_UNKNOWN_STRING, parameter::stack_type = "T"

  COOP_REAL, parameter::fwhm_arcmin = 30.d0
  COOP_REAL, parameter::fwhm_in = 20.d0
  COOP_UNKNOWN_STRING, parameter::prefix = "hslr5f30n1024/"
  COOP_UNKNOWN_STRING, parameter::postfix =   "_010a_1024.fits"

  COOP_UNKNOWN_STRING, parameter::mapdir = "/mnt/scratch-lustre/zqhuang/scratch-3month/zqhuang/"
  COOP_REAL,parameter::fwhm = coop_SI_arcmin * sqrt(fwhm_arcmin**2-fwhm_in**2)
  COOP_REAL, parameter::threshold = 0
  COOP_INT, parameter::mmax = 0
  COOP_REAL, parameter::dr = coop_SI_arcmin * max(fwhm_arcmin/2.d0, 10.d0)
  COOP_INT, parameter::n = nint(5.d0*coop_SI_degree/dr)

  COOP_UNKNOWN_STRING, parameter::imap_file  = "planck14/dx11_v2_smica_int_cmb"//postfix
  COOP_UNKNOWN_STRING, parameter::polmap_file  = "planck14/dx11_v2_smica_pol_case3_cmb"//postfix
  COOP_UNKNOWN_STRING, parameter::imask_file  = "planck14/dx11_v2_smica_int_mask"//postfix
  COOP_UNKNOWN_STRING, parameter::polmask_file  ="planck14/dx11_v2_smica_pol_mask"//postfix


  type(coop_healpix_maps)::polmask, imask, noise, imap, polmap, tmpmap
  type(coop_healpix_patch)::patch_s, patch_n
  integer,parameter::scan_nside = 4
  integer run_id, i, nlines, n_larger_chisq, ind
  COOP_REAL diff(0:n), chisq, prob, hdir(2), kdata, bdata, kmean, bmean, cov(2,2)
  COOP_STRING::fmt, fr_file, log_file, fig_file
  type(coop_list_integer)::listpix
  type(coop_list_real)::listangle
  type(coop_file) fp
  type(coop_asy) fig
  COOP_REAL,dimension(n_sim)::ksim, bsim
  COOP_REAL junk(6)
  call coop_MPI_init()
  if(iargc() .ge. 1)then
     run_id = coop_str2int(coop_InputArgs(1))
  else
     run_id = coop_MPI_Rank()
  endif

  if(run_id .ge.  scan_nside**2*12)then
     write(*,*) "run id must not exceed ", scan_nside**2*12 - 1
     call coop_MPI_Abort()
  endif
  call pix2ang_ring(scan_nside, run_id, hdir(1), hdir(2))

  !!read mask and map
  call imask%read(imask_file, nmaps_wanted = 1, spin = (/ 0 /) )
!  call polmask%read(polmask_file, nmaps_wanted = 1)

  call patch_n%init(stack_type, n, dr, mmax = mmax)
  patch_s = patch_n

  fr_file = prefix//stack_type//"_on_"//spot_type//"_fr_"//COOP_STR_OF(scan_nside)//"_"//COOP_STR_OF(run_id)//".txt"
  fig_file = prefix//stack_type//"_on_"//spot_type//"_fig_"//COOP_STR_OF(scan_nside)//"_"//COOP_STR_OF(run_id)//".txt"
  log_file = prefix//stack_type//"_on_"//spot_type//"_log_"//COOP_STR_OF(scan_nside)//"_"//COOP_STR_OF(run_id)//".txt"

  fmt = "("//trim(coop_num2str(n+1))//"G16.7)"
  ind = 0
  if(coop_file_exists(trim(fr_file)))then
     nlines = coop_file_numlines(trim(fr_file))
     if(nlines .gt. 0)then
        call fp%open(trim(fr_file), "r")
        do i = 1, nlines
           read(fp%unit, trim(fmt)) diff
           ind = ind + 1
           call coop_linear_least_square_fit(n+1, patch_s%r, diff, ksim(ind), bsim(ind))
           if(ind .ge. n_sim) exit
        enddo
        call fp%close()
        write(*,*) "Loaded "//trim(coop_num2str(ind))//" lines from checkpoint"
     endif
  endif


  call fp%open(trim(fr_file), "a")
  do while(ind .lt. n_sim)
     ind = ind + 1
     call load_imap(ind)
     noise%map = imap%map
     call noise%get_listpix(listpix, listangle, spot_type, threshold, imask)

     call imap%stack_north_south(patch_n, patch_s, listpix, listangle, hdir, imask)
     call patch_n%get_radial_profile(imap = 1, m = 0)
     call patch_s%get_radial_profile(imap = 1, m = 0)
     diff = patch_n%fr(:, 0, 1) - patch_s%fr(:, 0, 1)
     call coop_linear_least_square_fit(n+1, patch_s%r, diff, ksim(ind), bsim(ind))
     write(fp%unit, trim(fmt)) diff
     flush(fp%unit)
  enddo
  call fp%close()

  if(coop_file_exists(trim(log_file)))then
     call fp%open(trim(log_file), "r")
     read(fp%unit, *) junk
     read(fp%unit, *) diff
     read(fp%unit, *) kdata, bdata
     call fp%close()
  else
     call imap%read(imap_file, nmaps_wanted = 1, spin = (/ 0 /) )
     imap%map(:,1) = imap%map(:,1)*imask%map(:,1)
    if(fwhm.ge.coop_SI_arcmin)    call imap%smooth(fwhm)
     call imap%get_listpix(listpix, listangle, spot_type, threshold, imask)
     call imap%convert2ring()
     call imap%stack_north_south(patch_n, patch_s, listpix, listangle, hdir, imask )
     call patch_s%get_radial_profile(imap = 1, m = 0)
     call patch_n%get_radial_profile(imap = 1, m = 0)
     diff = patch_n%fr(:, 0, 1) - patch_s%fr(:, 0, 1)
     call coop_linear_least_square_fit(n+1, patch_s%r, diff, kdata, bdata)
  endif
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
  write(fp%unit, fmt) diff
  write(fp%unit, "(2E16.7)") kdata, bdata
  call fp%close()
     

  ksim = ksim*1.e6
  bsim = bsim*1.e6
  kdata = kdata*1.e6
  bdata = bdata*1.e6

  call fig%open(trim(fig_file))
  call fig%init(ylabel = "$\delta T (\mu K) $", xlabel = "$dT /dr (\mu K/{\rm rad})$")
  call coop_asy_dots(fig, ksim, bsim, "black")
  call coop_asy_dot(fig, kdata, bdata, "red", "x")
  call coop_asy_label(fig, "x : data",  fig%xmin+(fig%xmax - fig%xmin)*0.12, fig%ymin + 0.86*(fig%ymax - fig%ymin), color = "red")
  call coop_asy_label(fig, "$\bullet$ : sims",  fig%xmin+(fig%xmax - fig%xmin)*0.15, fig%ymin + 0.82*(fig%ymax - fig%ymin), color = "black")
  call coop_asy_label(fig, "direction: $l ="//trim(coop_num2str(nint(hdir(2)/coop_SI_degree)))//"^{\circ},  b="//trim(coop_num2str(nint((coop_pio2-hdir(1))/coop_SI_degree)))//"^{\circ}$",  fig%xmin+(fig%xmax - fig%xmin)*0.2, fig%ymin + 0.78*(fig%ymax - fig%ymin), color = "black")
  call fig%close()
  call coop_MPI_Finalize()
contains

  subroutine load_imap(i)
    COOP_INT i
    call imap%read(trim(sim_file_name_cmb_imap(i)), spin = (/ 0 /), nmaps_wanted = 1  )
    call noise%read(trim(sim_file_name_noise_imap(i)), spin = (/ 0 /) , nmaps_wanted = 1 )
    imap%map(:, 1) = (imap%map(:, 1) + noise%map(:, 1))*imask%map(:, 1)
    if(fwhm.ge.coop_SI_arcmin)    call imap%smooth(fwhm)
  end subroutine load_imap


  subroutine load_polmap(i)
    COOP_INT i
    call polmap%read(trim(sim_file_name_cmb_polmap(i)), spin = (/2 , 2 /) , nmaps_wanted = 2  )
    call noise%read(trim(sim_file_name_noise_polmap(i)), spin = (/2 , 2 /) , nmaps_wanted = 2 )
    polmap%map(:, 1) = (polmap%map(:, 1) + noise%map(:, 1))*polmask%map(:, 1)
    polmap%map(:, 2) = (polmap%map(:, 2) + noise%map(:, 2))*polmask%map(:, 1)
    if(fwhm.ge.coop_SI_arcmin)call polmap%smooth(fwhm)
  end subroutine load_polmap


  function sim_file_name_cmb_imap(i)
    COOP_INT i
    COOP_STRING sim_file_name_cmb_imap
    sim_file_name_cmb_imap = mapdir//"cmb/int/dx11_v2_smica_int_cmb_mc_"//trim(coop_Ndigits(i-1, 5))//postfix
  end function sim_file_name_cmb_imap

  function sim_file_name_noise_imap(i)
    COOP_INT i
    COOP_STRING sim_file_name_noise_imap
    sim_file_name_noise_imap = mapdir//"noise/int/dx11_v2_smica_int_noise_mc_"//trim(coop_Ndigits(i-1, 5))//postfix
  end function sim_file_name_noise_imap


  function sim_file_name_cmb_polmap(i)
    COOP_INT i
    COOP_STRING sim_file_name_cmb_polmap
    sim_file_name_cmb_polmap = mapdir//"cmb/pol/dx11_v2_smica_pol_case3_cmb_mc_"//trim(coop_Ndigits(i-1, 5))//postfix
  end function sim_file_name_cmb_polmap

  function sim_file_name_noise_polmap(i)
    COOP_INT i
    COOP_STRING sim_file_name_noise_polmap
    sim_file_name_noise_polmap = mapdir//"noise/pol/dx11_v2_smica_pol_case3_noise_mc_"//trim(coop_Ndigits(i-1, 5))//postfix
  end function sim_file_name_noise_polmap


end program test
