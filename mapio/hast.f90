program hastack_prog
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
  COOP_REAL, parameter::patch_size = 4.d0*coop_SI_degree

  

  COOP_UNKNOWN_STRING, parameter::prefix = "ha_r4f30n1024/"
    COOP_INT, parameter::mmax = 4
  COOP_REAL, parameter::fwhm_arcmin = 30.d0
  COOP_REAL, parameter::fwhm_in = 10.d0
  COOP_UNKNOWN_STRING, parameter::postfix =   "_010a_1024.fits"

  COOP_UNKNOWN_STRING, parameter::mapdir = "/mnt/scratch-lustre/zqhuang/scratch-3month/zqhuang/"
  COOP_REAL,parameter::fwhm = coop_SI_arcmin * sqrt(fwhm_arcmin**2-fwhm_in**2)
  COOP_REAL, parameter::threshold = 0
  COOP_REAL, parameter::dr = coop_SI_arcmin * max(fwhm_arcmin/4.d0, fwhm_in)
  COOP_INT, parameter::n = nint(patch_size/dr)

  COOP_UNKNOWN_STRING, parameter::imap_file  = "planck14/dx11_v2_smica_int_cmb"//postfix
  COOP_UNKNOWN_STRING, parameter::polmap_file  = "planck14/dx11_v2_smica_pol_case3_cmb"//postfix
  COOP_UNKNOWN_STRING, parameter::imask_file  = "planck14/dx11_v2_smica_int_mask"//postfix
  COOP_UNKNOWN_STRING, parameter::polmask_file  ="planck14/dx11_v2_smica_pol_mask"//postfix

  type(coop_healpix_maps)::polmask, imask, noise, imap, polmap, tmpmap
  type(coop_healpix_patch)::patch_s, patch_n
  integer,parameter::scan_nside = 4
    integer,parameter::scan_npix = scan_nside**2*6
  integer run_id, i, ind, j
  COOP_REAL   hdir(2)
  COOP_STRING::fr_file
  type(coop_list_integer)::listpix
  type(coop_list_real)::listangle
  type(coop_file) fp
  COOP_REAL,dimension(n_sim)::ksim, bsim
  COOP_REAL junk(6)
  call coop_MPI_init()
  if(iargc() .ge. 1)then
     run_id = coop_str2int(coop_InputArgs(1))
  else
     run_id = coop_MPI_Rank()
  endif

  if(run_id .ge.  scan_npix)then
     write(*,*) "run id must not exceed ", scan_nside**2*12 - 1
     call coop_MPI_Abort()
  endif
  call pix2ang_ring(scan_nside, run_id, hdir(1), hdir(2))

  !!read masks
  call imask%read(imask_file, nmaps_wanted = 1, spin = (/ 0 /) )
  if(stack_type .ne. "T")then
     call polmask%read(polmask_file, nmaps_wanted = 1, spin = (/ 0 /) )
  endif

  
  call patch_n%init(stack_type, n, dr, mmax = mmax)
  patch_s = patch_n

  fr_file = prefix//stack_type//"_on_"//spot_type//"_fr_"//COOP_STR_OF(scan_nside)//"_"//COOP_STR_OF(run_id)//".txt"


  ind = -1
  if(coop_file_exists(trim(fr_file)))then
     call fp%open(trim(fr_file), "ru")
     do
        read(fp%unit, ERR=100, END=100) i
        read(fp%unit, ERR=100, END=100) patch_n%fr
        read(fp%unit, ERR=100, END=100) patch_s%fr
        if(i.ne.ind+1) call cooP_MPI_Abort("fr file broken")
        ind = i
        if(ind .ge. n_sim) exit
     enddo
     call fp%close()
     write(*,*) "Loaded "//trim(coop_num2str(ind+1))//" lines from checkpoint"
  endif
100 call fp%open(trim(fr_file), "u")
  do i=0, ind
     read(fp%unit, ERR=100, END=100) j
     read(fp%unit, ERR=100, END=100) patch_n%fr
     read(fp%unit, ERR=100, END=100) patch_s%fr
     if(j.ne.i)call cooP_MPI_Abort("fr file broken")
  enddo
  do while(ind .lt. n_sim)
     ind = ind + 1
     call load_imap(ind)
     select case(stack_type)
     case("T")
        noise%map = imap%map
        call noise%get_listpix(listpix, listangle, spot_type, threshold, imask)
        call imap%stack_north_south(patch_n, patch_s, listpix, listangle, hdir, imask)
     case default
        call imap%get_listpix(listpix, listangle, spot_type, threshold, imask)
        call load_polmap(ind)
        call polmap%stack_north_south(patch_n, patch_s, listpix, listangle, hdir, polmask)
     end select

     call patch_n%get_all_radial_profiles()
     call patch_s%get_all_radial_profiles()
     write(fp%unit) ind
     write(fp%unit) patch_n%fr
     write(fp%unit) patch_s%fr
     flush(fp%unit)
  enddo
  call fp%close()
  call coop_MPI_Finalize()

contains

  subroutine load_imap(i)
    COOP_INT i, nm
    COOP_INT,dimension(:),allocatable::spin
    if(spot_type .eq. "PTmax" .or. spot_type .eq. "Tmax_QTUTOrient")then
       nm = 3
    else
       nm = 1
    endif
    allocate(spin(nm))
    spin(1) = 0
    if(nm.gt.1) spin(2:3) = 2
    if(i.eq.0)then
       call imap%read(filename = trim(imap_file), nmaps_wanted = nm , spin = spin , nmaps_to_read = 1 )
       imap%map(:, 1) = imap%map(:, 1)*imask%map(:, 1)
       if(fwhm.ge.coop_SI_arcmin)    call imap%smooth(fwhm)
       noise = imap
    else
       call imap%read(trim(sim_file_name_cmb_imap(i)), nmaps_wanted = nm , spin = spin , nmaps_to_read = 1 )
       call noise%read(trim(sim_file_name_noise_imap(i)), nmaps_wanted = nm , spin = spin, nmaps_to_read = 1 )
       imap%map(:, 1) = (imap%map(:, 1) + noise%map(:, 1))*imask%map(:, 1)
       if(fwhm.ge.coop_SI_arcmin)    call imap%smooth(fwhm)
    endif
    deallocate(spin)
    if(nm.gt.1)call imap%iqu2TQTUT()
  end subroutine load_imap


  subroutine load_polmap(i)
    COOP_INT i
    if(i.eq.0)then
       call polmap%read(trim(polmap_file), spin = (/2 , 2 /) , nmaps_wanted = 2  )
       polmap%map(:, 1) = polmap%map(:, 1)*polmask%map(:, 1)
       polmap%map(:, 2) = polmap%map(:, 2)*polmask%map(:, 1)
       if(fwhm.ge.coop_SI_arcmin)call polmap%smooth(fwhm)
    else
       call polmap%read(trim(sim_file_name_cmb_polmap(i)), spin = (/2 , 2 /) , nmaps_wanted = 2  )
       call noise%read(trim(sim_file_name_noise_polmap(i)), spin = (/2 , 2 /) , nmaps_wanted = 2 )
       polmap%map(:, 1) = (polmap%map(:, 1) + noise%map(:, 1))*polmask%map(:, 1)
       polmap%map(:, 2) = (polmap%map(:, 2) + noise%map(:, 2))*polmask%map(:, 1)
       if(fwhm.ge.coop_SI_arcmin)call polmap%smooth(fwhm)
    endif
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

end program hastack_prog
