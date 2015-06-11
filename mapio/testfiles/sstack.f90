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
  logical,parameter::do_calibration = .false.  
  COOP_INT, parameter::n_sim = 1000
  COOP_UNKNOWN_STRING, parameter::color_table = "Rainbow"
  COOP_SHORT_STRING::spot_type, stack_type, threshold_input
  
  COOP_REAL, parameter::patch_size = 2.d0*sin(coop_SI_degree)
  COOP_INT, parameter::n = 36
  COOP_REAL, parameter::dr = patch_size/n
  COOP_REAL::threshold = 1.  

  COOP_UNKNOWN_STRING, parameter::cs_method = "smica"
  COOP_UNKNOWN_STRING, parameter::pol_case = "case1"
  
  COOP_UNKNOWN_STRING, parameter::prefix = "sr/"//cs_method//"_"//pol_case//"_"
  COOP_INT, parameter::mmax = 4
  COOP_REAL, parameter::fwhm_arcmin = 15.d0
  COOP_REAL, parameter::fwhm_in = 10.d0
  COOP_UNKNOWN_STRING, parameter::input_resolution =   "_010a_1024"
  COOP_UNKNOWN_STRING, parameter::postfix =  input_resolution//".fits"
  COOP_UNKNOWN_STRING, parameter::polpost = "_hp_20_40"//postfix

  COOP_STRING::allprefix
  COOP_UNKNOWN_STRING, parameter::mapdir = "/mnt/scratch-lustre/zqhuang/scratch-3month/zqhuang/"
  COOP_REAL,parameter::fwhm = coop_SI_arcmin * sqrt(fwhm_arcmin**2-fwhm_in**2)
  
  COOP_UNKNOWN_STRING, parameter::imap_file  = "planck14/dx11_v2_"//cs_method//"_int_cmb"//postfix
  COOP_UNKNOWN_STRING, parameter::polmap_file  = "planck14/dx11_v2_"//cs_method//"_pol_"//pol_case//"_cmb"//polpost
  COOP_UNKNOWN_STRING, parameter::imask_file  = "planck14/dx11_v2_common_int_mask"//postfix
  COOP_UNKNOWN_STRING, parameter::polmask_file  ="planck14/dx11_v2_common_pol_mask"//postfix

  type(coop_healpix_maps)::polmask, imask, inoise, polnoise, imap, polmap, tmpmap
  type(coop_healpix_patch)::patch
  COOP_INT i, j, ind, nmaps_temp, l, lmax
  COOP_STRING::fr_file
  COOP_REAL::cl13, cl14
  type(coop_list_integer)::listpix
  type(coop_list_real)::listangle
  type(coop_file) fp
  logical::iload = .false.
  logical::polload = .false.
  COOP_SINGLE,dimension(:),allocatable::window
  logical::iload_first_call  = .true.
  logical::polload_first_call = .true.

  call coop_MPI_init()
  spot_type = trim(coop_InputArgs(1))
  stack_type = trim(coop_InputArgs(2))
  threshold_input = trim(coop_InputArgs(3))
  if(trim(threshold_input).ne."")then
     read(threshold_input,*)threshold
  endif
  if(trim(spot_type) .eq. "" .or. trim(stack_type).eq."")then
     print*, "Syntax:"
     print*, "./SST Tmax  T [nu]"
     print*, "./SST Tmax_QTUTOrient  T [nu]"     
     print*, "./SST Tmax  QrUr [nu]"
     print*, "./SST Tmax  QU [nu]"
     print*, "./SST Tmax_QTUTOrient QU [nu]"
     print*, "./SST PTmax QU [nu]"
     print*, "./SST Pmax QU [nu]"
     stop
  endif

  lmax  = 2000
  write(*,*) "Using lmax = "//COOP_STR_OF(lmax)
  allocate(window(0:lmax))
  window(0:1) = 0.d0  
  if(do_calibration)then
     call fp%open("calcls.txt", 'r')
     do l = 2, lmax
        read(fp%unit, *) i, cl13, cl14
        if(i.ne.l) stop "calcls.txt file broken"
        window(l) = sqrt(cl14/cl13)
     enddo
     call fp%close()
  else
     window(2:lmax) = 1.d0
  endif
  
  call imask%read(imask_file, nmaps_wanted = 1, spin = (/ 0 /) )  
  call polmask%read(polmask_file, nmaps_wanted = 1, spin = (/ 0 /) )
  
  call patch%init(trim(stack_type), n, dr, mmax = mmax)
  if(threshold .gt. -2.d0)then
     allprefix = prefix//trim(stack_type)//"_on_"//trim(spot_type)//"_fr_"//COOP_STR_OF(nint(patch_size/coop_SI_degree))//"deg"//input_resolution//"_smooth"//COOP_STR_OF(nint(fwhm_arcmin))//"_nu"//COOP_STR_OF(nint(threshold))     
  else
     allprefix = prefix//trim(stack_type)//"_on_"//trim(spot_type)//"_fr_"//COOP_STR_OF(nint(patch_size/coop_SI_degree))//"deg"//input_resolution//"_smooth"//COOP_STR_OF(nint(fwhm_arcmin))
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
     iload = .false.
     polload = .false.
     ind = ind + 1
     if(mod(ind, 5).eq.0)write(*,*) "Stacking map #"//COOP_STR_OF(ind)
     select case(trim(spot_type))
     case("Tmax", "Tmax_QTUTOrient", "PTmax")
        call load_imap(ind)
        call imap%get_listpix(listpix, listangle, trim(spot_type), threshold, imask)
     case("Pmax")
        call load_polmap(ind)
        call polmap%get_listpix(listpix, listangle, trim(spot_type), threshold, polmask)
     case("default")
        print*, trim(spot_type)
        stop "Not supported"
     end select
     select case(trim(stack_type))
     case("T")
        call load_imap(ind)
        call imap%stack_with_listpix(patch, listpix, listangle, imask)
     case("QU", "QrUr")
        call load_polmap(ind)
        call polmap%stack_with_listpix(patch, listpix, listangle, polmask)
     case default
        print*, trim(stack_type)
        stop "Not supported"
     end select
     call patch%get_all_radial_profiles()
     write(fp%unit) ind
     write(fp%unit) patch%fr
     if(mod(ind, 30).eq. 0) flush(fp%unit)
  enddo
  call fp%close()
  call coop_MPI_Finalize()

contains

  subroutine load_imap(i)
    COOP_INT i, nm
    COOP_INT,dimension(:),allocatable::spin
    if(iload)return
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
    else
       if(iload_first_call)then
          call imap%read(trim(sim_file_name_cmb_imap(i)), nmaps_wanted = nm , spin = spin , nmaps_to_read = 1 )
          call inoise%read(trim(sim_file_name_noise_imap(i)), nmaps_wanted = nm , spin = spin, nmaps_to_read = 1 )
          iload_first_call = .false.
       else
          call imap%read(trim(sim_file_name_cmb_imap(i)),  nmaps_to_read = 1, known_size = .true.)
          call inoise%read(trim(sim_file_name_noise_imap(i)), nmaps_to_read = 1, known_size = .true.)          
       endif
       imap%map(:, 1) = (imap%map(:, 1) + inoise%map(:, 1))*imask%map(:, 1)
    endif
    deallocate(spin)
    call do_smooth_imap(imap, i)        
    if(nm.gt.1)call imap%iqu2TQTUT( idone = (fwhm .gt. coop_SI_arcmin .or. do_calibration)  )
  end subroutine load_imap


  subroutine load_polmap(i)
    COOP_INT i
    if(polload)return
    if(i.eq.0)then
       call polmap%read(trim(polmap_file), spin = (/2 , 2 /) , nmaps_wanted = 2  )
       polmap%map(:, 1) = polmap%map(:, 1)*polmask%map(:, 1)
       polmap%map(:, 2) = polmap%map(:, 2)*polmask%map(:, 1)
    else
       if(polload_first_call)then
          call polmap%read(trim(sim_file_name_cmb_polmap(i)), spin = (/2 , 2 /) , nmaps_wanted = 2  )
          call polnoise%read(trim(sim_file_name_noise_polmap(i)), spin = (/2 , 2 /) , nmaps_wanted = 2 )
          polload_first_call = .false.
       else
          call polmap%read(trim(sim_file_name_cmb_polmap(i)), nmaps_to_read = 2 , known_size = .true.)
          call polnoise%read(trim(sim_file_name_noise_polmap(i)), nmaps_to_read = 2, known_size = .true.)
       endif       
       polmap%map(:, 1) = (polmap%map(:, 1) + polnoise%map(:, 1))*polmask%map(:, 1)
       polmap%map(:, 2) = (polmap%map(:, 2) + polnoise%map(:, 2))*polmask%map(:, 1)
    endif
    call do_smooth_polmap(polmap, i)    
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
    sim_file_name_cmb_polmap = mapdir//"cmb/pol/dx11_v2_"//cs_method//"_pol_"//pol_case//"_cmb_mc_"//trim(coop_Ndigits(i-1, 5))//polpost
  end function sim_file_name_cmb_polmap

  function sim_file_name_noise_polmap(i)
    COOP_INT i
    COOP_STRING sim_file_name_noise_polmap
    sim_file_name_noise_polmap = mapdir//"noise/pol/dx11_v2_"//cs_method//"_pol_"//pol_case//"_noise_mc_"//trim(coop_Ndigits(i-1, 5))//polpost
  end function sim_file_name_noise_polmap


  subroutine do_smooth_imap(mymap, i)
    type(coop_healpix_maps)mymap
    integer i
    if(fwhm .gt. coop_SI_arcmin .or. do_calibration)then
       if(do_calibration .and. i.gt.0)then
          call mymap%smooth_with_window(fwhm = fwhm, lmax = lmax, window = window, index_list = (/ 1 /) )          
       else
          call mymap%smooth(fwhm = fwhm, l_upper = lmax, l_lower = 2, index_list = (/ 1 /))
       endif
    endif
  end subroutine do_smooth_imap
  

  subroutine do_smooth_polmap(mymap, i)
    type(coop_healpix_maps)mymap
    integer i
    if(fwhm .gt. coop_SI_arcmin) call mymap%smooth(fwhm = fwhm, index_list = (/ 1 , 2 /) , l_upper = lmax, l_lower = 2)
  end subroutine do_smooth_polmap

end program hastack_prog
