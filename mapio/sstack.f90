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
  COOP_INT, parameter::n_sim = 100
  COOP_UNKNOWN_STRING, parameter::color_table = "Rainbow"
  COOP_SHORT_STRING::spot_type, stack_type
  COOP_REAL, parameter::patch_size = 2.d0*coop_SI_degree
  COOP_UNKNOWN_STRING, parameter::cs_method = "commander"
  COOP_UNKNOWN_STRING, parameter::pol_case = "case5"
  
  COOP_UNKNOWN_STRING, parameter::prefix = "sst_"//cs_method//"_"//pol_case//"/"
  COOP_INT, parameter::mmax = 4
  COOP_REAL, parameter::fwhm_arcmin = 30.d0
  COOP_REAL, parameter::fwhm_in = 10.d0
  COOP_UNKNOWN_STRING, parameter::input_resolution =   "_010a_1024"
  COOP_UNKNOWN_STRING, parameter::postfix =  input_resolution//".fits"

  COOP_STRING::allprefix
  COOP_UNKNOWN_STRING, parameter::mapdir = "/mnt/scratch-lustre/zqhuang/scratch-3month/zqhuang/"
  COOP_REAL,parameter::fwhm = coop_SI_arcmin * sqrt(fwhm_arcmin**2-fwhm_in**2)
  COOP_REAL, parameter::threshold = 0
  COOP_REAL, parameter::dr = coop_SI_arcmin * max(fwhm_arcmin/4.d0, 5.d0)
  COOP_INT, parameter::n = nint(patch_size/dr)

  COOP_UNKNOWN_STRING, parameter::imap_file  = "planck14/dx11_v2_"//cs_method//"_int_cmb"//postfix
  COOP_UNKNOWN_STRING, parameter::polmap_file  = "planck14/dx11_v2_"//cs_method//"_pol_"//pol_case//"_cmb"//postfix
  COOP_UNKNOWN_STRING, parameter::imask_file  = "planck14/dx11_v2_common_int_mask"//postfix
  COOP_UNKNOWN_STRING, parameter::polmask_file  ="planck14/dx11_v2_common_pol_mask"//postfix

  type(coop_healpix_maps)::polmask, imask, noise, imap, polmap, tmpmap
  type(coop_healpix_patch)::patch
  COOP_INT i, j, ind, nmaps_temp, l, lmax
  COOP_STRING::fr_file
  type(coop_list_integer)::listpix
  type(coop_list_real)::listangle
  type(coop_file) fp
  COOP_SINGLE,dimension(:),allocatable::window
  logical highpass
  call coop_MPI_init()
  spot_type = trim(coop_InputArgs(1))
  stack_type = trim(coop_InputArgs(2))
  if(trim(spot_type) .eq. "" .or. trim(stack_type).eq."")then
     print*, "Syntax:"
     print*, "./SST Tmax  T"
     print*, "./SST Tmax  QrUr [highpass]"
     print*, "./SST Tmax  QU [highpass]"     
     print*, "./SST Tmax_QTUTOrient QU [highpass]"
     print*, "./SST PTmax QU [highpass]"
     print*, "where [highpass] can be T or F or ignored (default F, pol only)"
     stop
  endif

  lmax  = ceiling(3.d0/(fwhm_arcmin*coop_SI_arcmin * coop_sigma_by_fwhm))
  write(*,*) "Using lmax = "//COOP_STR_OF(lmax)
  highpass = (trim(coop_InputArgs(3)) .eq. "T")
  call imask%read(imask_file, nmaps_wanted = 1, spin = (/ 0 /) )
  if(highpass)then
     allocate(window(0:lmax))
     do l= 0, lmax
        window(l) = high_pass_window(l, 20, 40) !!default 20, 40
     enddo
  endif
  
  if(trim(stack_type) .ne. "T")then
     call polmask%read(polmask_file, nmaps_wanted = 1, spin = (/ 0 /) )
  endif
  
  call patch%init(trim(stack_type), n, dr, mmax = mmax)
  if(fwhm .gt. coop_SI_arcmin)then
     allprefix = prefix//trim(stack_type)//"_on_"//trim(spot_type)//"_fr_"//COOP_STR_OF(nint(patch_size/coop_SI_degree))//"deg"//input_resolution//"_smooth"//COOP_STR_OF(nint(fwhm_arcmin))
  else
     allprefix = prefix//trim(stack_type)//"_on_"//trim(spot_type)//"_fr_"//COOP_STR_OF(nint(patch_size/coop_SI_degree))//"deg"//input_resolution//"_nosmooth"
  endif
  if(highpass)then
     allprefix = trim(allprefix)//"_highpass"
  endif
  call fp%open(trim(allprefix)//"_info.txt", "w")
  write(fp%unit,*) n, patch%nmaps, dr/coop_SI_arcmin
  call fp%close()

  fr_file = trim(allprefix)//".dat"
  ind = -1
  if(.not. coop_file_exists(trim(fr_file)))goto 200
  call fp%open(trim(fr_file), "ru")
  do
     read(fp%unit, ERR=100, END=100) i
     read(fp%unit, ERR=100, END=100) patch%fr
     if(i.ne.ind+1) call cooP_MPI_Abort("fr file broken")
     ind = i
     if(ind .ge. n_sim) exit
  enddo
100 write(*,*) "Loaded "//trim(coop_num2str(ind+1))//" maps from checkpoint"
  call fp%close()
200 call fp%open(trim(fr_file), "u")
  do i=0, ind
     read(fp%unit, ERR=100, END=100) j
     read(fp%unit, ERR=100, END=100) patch%fr
     if(j.ne.i)call cooP_MPI_Abort("fr file broken")
  enddo
  do while(ind .lt. n_sim)
     ind = ind + 1
     write(*,*) "Stacking map #"//COOP_STR_OF(ind)
     call load_imap(ind)
     select case(trim(stack_type))
     case("T")
        noise%map = imap%map
        call noise%get_listpix(listpix, listangle, trim(spot_type), threshold, imask)
        call imap%stack_with_listpix(patch, listpix, listangle, imask)
     case default
        call imap%get_listpix(listpix, listangle, trim(spot_type), threshold, imask)
        call load_polmap(ind)
        call polmap%stack_with_listpix(patch, listpix, listangle, polmask)
     end select
     call patch%get_all_radial_profiles()
     write(fp%unit) ind
     write(fp%unit) patch%fr
     flush(fp%unit)
  enddo
  call fp%close()
  call coop_MPI_Finalize()

contains

  subroutine load_imap(i)
    COOP_INT i, nm
    COOP_INT,dimension(:),allocatable::spin
    if(trim(spot_type) .eq. "PTmax" .or. trim(spot_type) .eq. "Tmax_QTUTOrient")then
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
       call do_smooth_map(imap)
       noise = imap
    else
       call imap%read(trim(sim_file_name_cmb_imap(i)), nmaps_wanted = nm , spin = spin , nmaps_to_read = 1 )
       call noise%read(trim(sim_file_name_noise_imap(i)), nmaps_wanted = nm , spin = spin, nmaps_to_read = 1 )
       imap%map(:, 1) = (imap%map(:, 1) + noise%map(:, 1))*imask%map(:, 1)
       call do_smooth_map(imap)
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
       call do_smooth_map(polmap)
    else
       call polmap%read(trim(sim_file_name_cmb_polmap(i)), spin = (/2 , 2 /) , nmaps_wanted = 2  )
       call noise%read(trim(sim_file_name_noise_polmap(i)), spin = (/2 , 2 /) , nmaps_wanted = 2 )
       polmap%map(:, 1) = (polmap%map(:, 1) + noise%map(:, 1))*polmask%map(:, 1)
       polmap%map(:, 2) = (polmap%map(:, 2) + noise%map(:, 2))*polmask%map(:, 1)
       call do_smooth_map(polmap)
    endif
  end subroutine load_polmap


  function sim_file_name_cmb_imap(i)
    COOP_INT i
    COOP_STRING sim_file_name_cmb_imap
    sim_file_name_cmb_imap = mapdir//"cmb/int/dx11_v2_"//cs_method//"_int_cmb_mc_"//trim(coop_Ndigits(i-1, 5))//postfix
  end function sim_file_name_cmb_imap

  function sim_file_name_noise_imap(i)
    COOP_INT i
    COOP_STRING sim_file_name_noise_imap
    sim_file_name_noise_imap = mapdir//"noise/int/dx11_v2_"//cs_method//"_int_noise_mc_"//trim(coop_Ndigits(i-1, 5))//postfix
  end function sim_file_name_noise_imap


  function sim_file_name_cmb_polmap(i)
    COOP_INT i
    COOP_STRING sim_file_name_cmb_polmap
    sim_file_name_cmb_polmap = mapdir//"cmb/pol/dx11_v2_"//cs_method//"_pol_"//pol_case//"_cmb_mc_"//trim(coop_Ndigits(i-1, 5))//postfix
  end function sim_file_name_cmb_polmap

  function sim_file_name_noise_polmap(i)
    COOP_INT i
    COOP_STRING sim_file_name_noise_polmap
    sim_file_name_noise_polmap = mapdir//"noise/pol/dx11_v2_"//cs_method//"_pol_"//pol_case//"_noise_mc_"//trim(coop_Ndigits(i-1, 5))//postfix
  end function sim_file_name_noise_polmap

  function high_pass_window(l, lstart, lend) result(w)
    COOP_SINGLE w
    COOP_INT l, lstart, lend
    if(l.le.lstart)then
       w = 0.
       return
    endif
    if(l.ge.lend)then
       w = 1.
       return
    endif
    w = sin(coop_pio2*dble(l - lstart)/dble(lend-lstart))**2
  end function high_pass_window

  subroutine do_smooth_map(mymap)
    type(coop_healpix_maps)mymap
    if(highpass .and. mymap%nmaps.ge.2 .and. all(mymap%spin.eq.2))then
       call mymap%smooth_with_window(fwhm = fwhm, lmax = lmax, window = window)          
    else
       if(fwhm .gt. coop_SI_arcmin) call mymap%smooth(fwhm = fwhm, l_upper = lmax, l_lower = 2)
    endif
  end subroutine do_smooth_map

end program hastack_prog
