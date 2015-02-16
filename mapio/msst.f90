program massive_stack
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX

  COOP_INT,parameter::num_nu = 5
  COOP_REAL, dimension(num_nu), parameter::nu = (/ 0.d0, 0.5d0, 1.d0, 1.5d0, 2.d0 /)
  COOP_INT, parameter::num_masks = 2
  COOP_STRING,dimension(num_masks)::sub_imask_file
  COOP_STRING::outputdir = "msst/"
  logical, parameter::do_nest = .true.

  COOP_STRING,parameter::peak_name = "$T$"
  COOP_INT,parameter::n = 36
  COOP_REAL,parameter::r_degree  = 2.d0
  COOP_REAL,parameter::dr = 2.d0*sin(r_degree*coop_SI_degree/2.d0)/n
  COOP_UNKNOWN_STRING,parameter::mapdir = "massffp8/"  
  
  COOP_STRING::cc_method 
  COOP_STRING::polcase
  COOP_STRING::stack_field_name
  COOP_INT::resol = 1024
  COOP_INT::fwhm
  COOP_STRING::output, line
  COOP_INT::n_sim
  COOP_STRING::orient_name   
  COOP_STRING::postfix
  COOP_STRING::imap_file, polmap_file, imask_file, polmask_file
  type(coop_stacking_options)::sto_max_all, sto_min_all  
  type(coop_stacking_options),dimension(num_masks)::sub_sto_max, sub_sto_min
  type(coop_stacking_options),dimension(num_nu, num_masks)::sto_max, sto_min
  type(coop_healpix_patch),dimension(num_nu, 0:num_masks)::patch_max, patch_min
  type(coop_healpix_maps)::imap, imask, polmask, inoise, polnoise, polmap
  type(coop_healpix_maps),dimension(num_masks)::sub_imask
  logical::iloaded = .false.
  logical::polloaded  = .false.
  type(coop_file)::fp
  COOP_INT i, j, in, im, nmaps, mmby2
  COOP_INT ind, ind_done
  COOP_REAL,dimension(:,:,:,:,:,:),allocatable::frall
  
  if(iargc() .lt. 5)then
     write(*,*) "Syntax:"
     write(*,*) "./MST cc_method resolution orient stack_field n_sim"
     write(*,*) "Examples:"     
     write(*,*) "./MST smica     1024 RANDOM  QU 1000"
     write(*,*) "./MST sevem      512 QTUT    T   500"
     write(*,*) "./MST nilc       512 QTUT    QU  500"
     write(*,*) "./MST commander 1024 RANDOM  T  1000"
     stop
  else
     output = trim(outputdir)//trim(coop_inputArgs(1))//"_"//trim(coop_inputArgs(2))//"_"//trim(coop_inputArgs(3))//"_"//trim(coop_inputArgs(4))//".dat"
  endif
  coop_healpix_warning = .false.
  cc_method = trim(coop_inputArgs(1))
  select case(trim(cc_method))
  case("smica","sevem","nilc")
     polcase= "pol_case1"
  case("commander")
     polcase = "pol_case5"
  case default
     write(*,*) "unknown cc_method = "//trim(cc_method)
     stop
  end select
  line = coop_inputArgs(2)
  read(line,*) resol  
  if(resol .ne. 1024 .and. resol .ne. 512) stop "only support resolution 512 and 1024"
  fwhm = 10240/resol
  
  select case(trim(coop_inputArgs(3)))
  case("QTUT")
     orient_name = "$(Q_T, U_T)$"
  case("RANDOM", "NULL")
     orient_name = "NULL"
  case default
     write(*,*) trim(coop_inputArgs(3))//" : unknown orientation"
     stop
  end select
  
  stack_field_name = trim(coop_inputArgs(4))
  if(trim(stack_field_name) .ne. "T" .and.   trim(stack_field_name) .ne. "QU")then
     stop "MST only support QU/T stacking"
  endif
  line = coop_inputArgs(5)
  read(line,*) n_sim

  postfix = "_"//trim(coop_ndigits(fwhm,3))//"a_"//trim(coop_ndigits(resol,4))//".fits"
  
  imap_file = "planck14/dx11_v2_"//trim(cc_method)//"_int_cmb"//trim(postfix)
  polmap_file = "planck14/dx11_v2_"//trim(cc_method)//"_"//trim(polcase)//"_cmb_hp_20_40"//trim(postfix)
  
  imask_file = "planck14/dx11_v2_common_int_mask"//trim(postfix)
  polmask_file = "planck14/dx11_v2_common_pol_mask"//trim(postfix)          

  sub_imask_file(1) = "planck14/HemAsym_north_int_mask"//trim(postfix)
  sub_imask_file(2) = "planck14/HemAsym_south_int_mask"//trim(postfix)


  call patch_max(1,0)%init(trim(stack_field_name), n, dr)
  nmaps = patch_max(1,0)%nmaps
  mmby2 = patch_max(1,0)%mmax/2
  
  allocate(frall(0:n, 0:mmby2, 1:nmaps, 1:num_nu, 0:num_masks,1:2)  )
  ind_done = -1  
  if(coop_file_exists(output))then
     call fp%open(output, "ru")
     do i = 0, n_sim
        read(fp%unit, ERR = 100, END = 100)  ind, frall(0:n, 0:mmby2, 1:nmaps, 1:num_nu, 0:num_masks,1:2)
        if(ind.ne.i) stop "data file broken"
        do j = 0, n
           write(*,*) j*dr, frall(j, :, 1, 1, 0, 1)
        enddo
        
        ind_done = ind
     enddo
100  call fp%close()
  endif

  coop_healpix_mask_tol = 0.5d0

  if(ind_done .ge. 0 .and. n_sim .gt. 0 )then
     print*, "loaded "//COOP_STR_OF(ind_done+1)//" stacked maps"
  endif

  if(ind_done .lt. n_sim)then
     do im = 0, num_masks     
        do in = 1, num_nu
           if(in.ne.1 .or. im.ne.0) patch_max(in, im) = patch_max(1, 0)
           patch_min(in, im) = patch_max(1, 0)        
        enddo
     enddo
     
     if(trim(orient_name).eq."NULL")then
        call sto_max_all%init(.true., peak_name, orient_name, nmaps = 1)
        call sto_min_all%init(.false., peak_name, orient_name, nmaps = 1)     
     else
        call sto_max_all%init(.true., peak_name, orient_name, nmaps = 3)  
        call sto_min_all%init(.false., peak_name, orient_name, nmaps = 3)     
     endif
     sto_max_all%threshold_option = 4
     sto_min_all%threshold_option = 4
     sto_max_all%nested = do_nest
     sto_min_all%nested = do_nest
     sto_max_all%I_lower_nu = nu(1)
     sto_min_all%I_upper_nu = -nu(1)
     do im = 1, num_masks
        sub_sto_max(im) = sto_max_all
        sub_sto_min(im) = sto_min_all
        do in = 1, num_nu
           sto_max(in,im) = sto_max_all
           sto_min(in,im) = sto_min_all
           sto_max(in,im)%I_lower_nu = nu(in)
           sto_min(in,im)%I_upper_nu = -nu(in)
           if(in.lt.num_nu)then
              sto_max(in, im)%I_upper_nu = nu(in+1)
              sto_min(in, im)%I_lower_nu = -nu(in+1)
           endif
        enddo
     enddo
     
     call imask%read(imask_file, nmaps_wanted = 1, spin = (/ 0 /) )
     imask%mask_npix = count(imask%map(:,1).gt. 0.5)
     do i=1, num_masks
        call sub_imask(i)%read(trim(sub_imask_file(i)), nmaps_wanted = 1, spin = (/ 0 /) )
     enddo
     if(do_nest)then
        call imask%convert2nested()
        do im=1, num_masks
           call sub_imask(im)%convert2nested()
        enddo
     endif
     if(trim(stack_field_name).eq."QU")then
        call polmask%read(polmask_file, nmaps_wanted = 1, spin = (/ 0 /) )
        if(do_nest) call polmask%convert2nested()
     endif
     call fp%open(output, "u")
     do ind = 1 , ind_done
        read(fp%unit)  i, frall(0:n, 0:mmby2, 1:nmaps, 1:num_nu, 0:num_masks,1:2)
     enddo
     do ind = ind_done+1, n_sim
        if(n_sim .gt. 0 )print*, "stacking map#"//COOP_STR_OF(ind)
        iloaded = .false.
        polloaded = .false.
        call find_peaks()
        call stack_map()
        call compute_fr()
        write(fp%unit) ind, frall(0:n, 0:mmby2, 1:nmaps, 1:num_nu, 0:num_masks,1:2)
     enddo
     call fp%close()     
  endif


500  continue
  
contains

  subroutine find_peaks()
    COOP_INT::im,in
    select case(sto_max_all%genre)
    case(coop_stacking_genre_Imax, coop_stacking_genre_Imax_Oriented)
       call load_imap(ind)
    case default
       stop "so far SST only support temperature peaks (oriented or nonoriented)"
    end select
    call imap%get_peaks(sto_max_all, mask = imask, restore = .not. do_nest)
    call imap%get_peaks(sto_min_all, mask = imask, restore = .not. do_nest)
    do im = 1, num_masks
       call sub_imask(im)%mask_peaks(sto_max_all, sub_sto_max(im))
       call sub_imask(im)%mask_peaks(sto_min_all, sub_sto_min(im))
       do in = 1, num_nu
          call sub_sto_max(im)%subset(sto_max(in, im))
          call sub_sto_min(im)%subset(sto_min(in, im))
       enddo
    enddo
  end subroutine find_peaks

  subroutine compute_fr()
    COOP_INT::in, im, imap

   
    do im = 1, num_masks
       do in=num_nu-1, 1, -1
          patch_max(in, im)%image = patch_max(in, im)%image + patch_max(in+1, im)%image
          patch_min(in, im)%image = patch_min(in, im)%image + patch_min(in+1, im)%image
          patch_max(in,im)%nstack = patch_max(in,im)%nstack + patch_max(in+1,im)%nstack
          patch_min(in,im)%nstack = patch_min(in,im)%nstack + patch_min(in+1,im)%nstack
       enddo
    enddo
    do in=1, num_nu
       patch_max(in, 0)%image = 0.d0
       patch_min(in, 0)%image = 0.d0
       patch_max(in, 0)%nstack = 0.d0
       patch_min(in, 0)%nstack = 0.d0
       do im = 1, num_masks
          patch_max(in, 0)%image = patch_max(in, 0)%image + patch_max(in, im)%image
          patch_min(in, 0)%image = patch_min(in, 0)%image + patch_min(in, im)%image
          patch_max(in, 0)%nstack = patch_max(in, 0)%nstack + patch_max(in, im)%nstack
          patch_min(in, 0)%nstack = patch_min(in, 0)%nstack + patch_min(in, im)%nstack          
          
       enddo
    enddo    
    do im = 0, num_masks
       do in = 1, num_nu
          do imap  = 1, nmaps
             patch_max(in, im)%image(:,:,imap) = patch_max(in, im)%image(:,:,imap) /max(patch_max(in, im)%nstack, 1.d0)
             patch_min(in, im)%image(:,:,imap) = patch_min(in, im)%image(:,:,imap) /max(patch_min(in, im)%nstack, 1.d0)             
          enddo
          call patch_max(in,im)%get_all_radial_profiles()
          call patch_min(in,im)%get_all_radial_profiles()
          frall(:,:,:,in,im,1) = patch_max(in,im)%fr*1.d6
          frall(:,:,:,in,im,2) = patch_min(in,im)%fr*1.d6    
       enddo
    enddo
  end subroutine compute_fr

  subroutine stack_map()
    COOP_INT in, im
    select case(trim(stack_field_name))
    case("T")
       call load_imap(ind)
       do im=1, num_masks
          do in=1, num_nu
             call imap%stack_on_peaks(sto_max(in, im), patch_max(in, im), mask = imask, norm = .false.)
             call imap%stack_on_peaks(sto_min(in, im), patch_min(in, im), mask = imask, norm = .false.)
          enddo
       enddo
    case("QU")
       call load_polmap(ind)
       do im=1, num_masks
          do in=1, num_nu
             call polmap%stack_on_peaks(sto_max(in, im), patch_max(in, im), mask = polmask, norm = .false.)
             call polmap%stack_on_peaks(sto_min(in, im), patch_min(in, im), mask = polmask, norm = .false.)
          enddo
       enddo
    case default
       stop "so far SST only support T and QU stacking"
    end select
  end subroutine stack_map

  subroutine load_imap(i)
    COOP_INT i
    if(iloaded) return
    if(i.eq.0)then
       call imap%read(imap_file, nmaps_wanted = sto_max_all%nmaps, nmaps_to_read = 1)
    else
       call imap%read(trim(mapdir)//"cmb/int/dx11_v2_"//trim(cc_method)//"_int_cmb_mc_"//trim(coop_Ndigits(i-1, 5))//trim(postfix), nmaps_to_read = 1, nmaps_wanted = sto_max_all%nmaps)
       call inoise%read(trim(mapdir)//"noise/int/dx11_v2_"//trim(cc_method)//"_int_noise_mc_"//trim(coop_Ndigits(i-1, 5))//trim(postfix), nmaps_to_read = 1, nmaps_wanted = sto_max_all%nmaps)
       imap%map(:, 1) = imap%map(:, 1) + inoise%map(:, 1)       
    endif
    if(imap%nmaps .eq. 3)then
       call imap%regularize_in_mask(imask_copy, 1)
       call imap%iqu2TQTUT()
    endif
    if(do_nest) call imap%convert2nested()
    iloaded = .true.
  end subroutine load_imap

  subroutine load_polmap(i)
    COOP_INT i
    if(polloaded) return
    if(i.eq.0)then
       call polmap%read(polmap_file, nmaps_wanted = 2, spin = (/ 2 , 2 /) )
    else
       call polmap%read(trim(mapdir)//"cmb/pol/dx11_v2_"//trim(cc_method)//"_"//trim(polcase)//"_cmb_mc_"//trim(coop_Ndigits(i-1, 5))//"_hp_20_40"//trim(postfix), nmaps_wanted = 2, spin = (/ 2, 2 /) )
       call polnoise%read(trim(mapdir)//"noise/pol/dx11_v2_"//trim(cc_method)//"_"//trim(polcase)//"_noise_mc_"//trim(coop_Ndigits(i-1, 5))//"_hp_20_40"//trim(postfix), nmaps_wanted = 2, spin = (/ 2 , 2 /) )
       polmap%map = polmap%map + polnoise%map
    endif
    if(do_nest) call polmap%convert2nested()
    polloaded = .true.
  end subroutine load_polmap

#else
  print*, "You need to install healpix"
#endif  
end program massive_stack


