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
  COOP_INT, parameter::n_sim = 30
  COOP_INT, parameter::mmax = 4  
  COOP_UNKNOWN_STRING, parameter::color_table = "Rainbow"
  COOP_SHORT_STRING::spot_type, stack_type, threshold_input

  
  COOP_UNKNOWN_STRING, parameter::prefix = "sr30/"
  COOP_REAL,parameter::r_degree = 2.d0
  COOP_REAL::threshold
  
  COOP_STRING::allprefix
  COOP_REAL,parameter::patch_size=2.d0*sin(r_degree*coop_SI_degree/2.d0)
  COOP_INT,parameter:: n_bins = 6
  COOP_INT, parameter:: n_per_bin = 8
  COOP_INT, parameter::n = n_bins * n_per_bin
  COOP_REAL, parameter::dr = patch_size/n

  logical, parameter::do_mask = .false.  
  COOP_UNKNOWN_STRING, parameter::imask_file  = "planck14/dx11_v2_common_int_mask_010a_1024.fits"
  COOP_UNKNOWN_STRING, parameter::polmask_file  ="planck14/dx11_v2_common_pol_mask_010a_1024.fits"

  type(coop_healpix_maps)::polmask, imask, imap, polmap, tmpmap
  type(coop_healpix_patch)::patch
  COOP_INT i, j, ind, nmaps_temp, l, lmax
  COOP_STRING::fr_file
  COOP_REAL::cl13, cl14
  type(coop_list_integer)::listpix
  type(coop_list_real)::listangle
  type(coop_file) fp
  logical::iload = .false.
  logical::polload = .false.
  logical::iload_first_call  = .true.
  logical::polload_first_call = .true.

  call coop_MPI_init()
  spot_type = trim(coop_InputArgs(1))
  stack_type = trim(coop_InputArgs(2))
  threshold_input = trim(coop_InputArgs(3))
  if(trim(threshold_input).ne."")then
     read(threshold_input,*)threshold
  else
     write(*,*) "No Threshold"
     threshold = 1.1e31
  endif
  if(trim(spot_type) .eq. "" .or. trim(stack_type).eq."")then
     print*, "Syntax:"
     print*, "./GetRad Tmax T [nu]"
     print*, "./GetRad Tmax_QTUTOrient T [nu]"     
     print*, "./GetRad Tmax QrUr [nu]"
     print*, "./GetRad Tmax QU [nu]"
     print*, "./GetRad Tmax_QTUTOrient QU [nu]"
     print*, "./GetRad PTmax QU [nu]"
     print*, "./GetRad Pmax QU [nu]"
     stop
  endif

  lmax  = 2000
  write(*,*) "Using lmax = "//COOP_STR_OF(lmax)
  if(do_mask)then
     call imask%read(imask_file, nmaps_wanted = 1, spin = (/ 0 /) )  
     call polmask%read(polmask_file, nmaps_wanted = 1, spin = (/ 0 /) )
  endif
  call patch%init(trim(stack_type), n, dr, mmax = mmax)
  if(abs(threshold).lt.6.d0)then
     allprefix = prefix//trim(stack_type)//"_on_"//trim(spot_type)//"_"//COOP_STR_OF(nint(r_degree))//"deg_nu"//COOP_STR_OF(nint(threshold))     
  else
     allprefix = prefix//trim(stack_type)//"_on_"//trim(spot_type)//"_"//COOP_STR_OF(nint(r_degree))//"deg_nuNone"
  endif

  call fp%open(trim(allprefix)//"_info.txt", "w")
  write(fp%unit,*) n, patch%nmaps, dr/coop_SI_arcmin, n_sim
  call fp%close()

  fr_file = trim(allprefix)//".dat"
  ind = 0
  if(.not. coop_file_exists(trim(fr_file)))goto 200
  call fp%open(trim(fr_file), "ru")
  do
     read(fp%unit, ERR=100, END=100) i
     read(fp%unit, ERR=100, END=100) patch%fr
     if(i.ne.ind+1) call cooP_MPI_Abort("fr file broken")
     ind = i
     if(ind .ge. n_sim) exit
  enddo
100 write(*,*) "Loaded "//trim(coop_num2str(ind))//" maps from checkpoint"
  call fp%close()
200 call fp%open(trim(fr_file), "u")
  do i=1, ind
     read(fp%unit, ERR=100, END=100) j
     read(fp%unit, ERR=100, END=100) patch%fr
     if(j.ne.i)call cooP_MPI_Abort("fr file broken")
  enddo
  do while(ind .lt. n_sim)
     iload = .false.
     polload = .false.
     ind = ind + 1
     write(*,*) "Stacking map #"//COOP_STR_OF(ind)
     select case(trim(spot_type))
     case("Tmax", "Tmax_QTUTOrient", "PTmax")
        call load_imap(ind)
        if(do_mask)then
           call imap%get_listpix(listpix, listangle, trim(spot_type), threshold, imask)
        else
           call imap%get_listpix(listpix, listangle, trim(spot_type), threshold)
        endif
     case("Pmax")
        call load_polmap(ind)
        if(do_mask)then
           call polmap%get_listpix(listpix, listangle, trim(spot_type), threshold, polmask)
        else
           call polmap%get_listpix(listpix, listangle, trim(spot_type), threshold)
        endif
     case("default")
        print*, trim(spot_type)
        stop "Not supported"
     end select
     select case(trim(stack_type))
     case("T")
        call load_imap(ind)
        if(do_mask)then
           call imap%stack_with_listpix(patch, listpix, listangle, imask)
        else
           call imap%stack_with_listpix(patch, listpix, listangle)           
        endif
     case("QU", "QrUr")
        call load_polmap(ind)
        if(do_mask)then
           call polmap%stack_with_listpix(patch, listpix, listangle, polmask)
        else
           call polmap%stack_with_listpix(patch, listpix, listangle)           
        endif
     case default
        print*, trim(stack_type)
        stop "Not supported"
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
    if(iload)return
    if(trim(spot_type) .eq. "PTmax" .or. trim(spot_type) .eq. "Tmax_QTUTOrient")then
       nm = 3
    else
       nm = 1
    endif
    allocate(spin(nm))
    spin(1) = 0
    if(nm.gt.1) spin(2:3) = 2
    if(iload_first_call)then
       call imap%read(trim(sim_file_name_cmb_imap(i)), nmaps_wanted = nm , spin = spin)
       iload_first_call = .false.
    else
       call imap%read(trim(sim_file_name_cmb_imap(i)),  nmaps_to_read = nm, known_size = .true.)
    endif
    deallocate(spin)
!    imap%map(:, 1) = imap%map(:, 1)*imask%map(:, 1)    
  end subroutine load_imap


  subroutine load_polmap(i)
    COOP_INT i
    if(polload)return
    if(polload_first_call)then
       call polmap%read(trim(sim_file_name_cmb_polmap(i)), spin = (/2 , 2 /) , nmaps_wanted = 2  )
       polload_first_call = .false.
    else
       call polmap%read(trim(sim_file_name_cmb_polmap(i)), nmaps_to_read = 2 , known_size = .true.)
    endif
!    polmap%map(:, 1) = polmap%map(:, 1)*polmask%map(:, 1)
!    polmap%map(:, 2) = polmap%map(:, 2)*polmask%map(:, 1)
  end subroutine load_polmap


  function sim_file_name_cmb_imap(i)
    COOP_INT i
    COOP_STRING sim_file_name_cmb_imap
    sim_file_name_cmb_imap ="massive/simu_TQTUT_"//trim(coop_ndigits(i-1, 5))//"_030a_n1024.fits"
  end function sim_file_name_cmb_imap


  function sim_file_name_cmb_polmap(i)
    COOP_INT i
    COOP_STRING sim_file_name_cmb_polmap
    sim_file_name_cmb_polmap = "massive/simu_pol_"//trim(coop_ndigits(i-1, 5))//"_030a_n1024.fits"
  end function sim_file_name_cmb_polmap

end program hastack_prog
