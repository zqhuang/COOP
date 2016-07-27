program coop_Xfaster
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use coop_fitsio_mod
  implicit none
#include "constants.h"
#define QB_IND(ib, icross, icl) (((ib-1)*numcross + icross - 1 )*numcls + icl)
#define CHIND(i)  trim(channel_name(i))
#define SIMIND(i) trim(sim_index_name(i))
#define CHSIMIND(ich, isim) CHIND(ich)//"_"//SIMIND(isim)
#define MUB_IND(ib, icl)  ((ib-1)*numcls + icl)

!!this is coop version of xfaster, written from scratch. 
  COOP_STRING::inifile, qb_output_root, cl_output_root, mub_output_root, Fisher_output_root, model_cl_file
  COOP_SHORT_STRING,dimension(:),allocatable::channel_name
  COOP_INT::sim_index_width, feedback
  COOP_STRING::action, map_genre, map_unit
  COOP_INT::ib, num_iterations, iter, ich
  COOP_REAL,parameter::fisher_threshold = 1.d-4
  COOP_SINGLE::map_maxval, map_minval
  type(coop_dictionary)::settings
  type(coop_file)::fp
  COOP_INT::num_ell_bins, num_channels, lmin_data, lmax_data, nmaps, numcls, numcross, num_qbs, matdim, num_mubs
  COOP_INT,dimension(:),allocatable::spin
  COOP_REAL,dimension(:),allocatable::qb, mub
  COOP_REAL junk
  COOP_STRING::fieldnames

  if(iargc().lt.1)then
     write(*,*) "Syntax:"
     write(*,*) "./XFASTER inifile"
     stop
  endif
  inifile = trim(coop_InputArgs(1))
  call coop_load_dictionary(inifile, settings)
  call coop_dictionary_lookup(settings, 'action', action)
  call coop_dictionary_lookup(settings, 'feedback', feedback, 1)
  call coop_dictionary_lookup(settings, 'map_genre', map_genre, 'IQU')
  call coop_dictionary_lookup(settings, 'qb_output_root', qb_output_root)
  call coop_dictionary_lookup(settings, 'mub_output_root', mub_output_root)
  call coop_dictionary_lookup(settings, 'fisher_output_root', fisher_output_root)

  call coop_dictionary_lookup(settings, 'model_cl_file', model_cl_file)
  call coop_dictionary_lookup(settings, 'cl_output_root', cl_output_root)
  call coop_dictionary_lookup(settings, 'num_iterations', num_iterations, 1)

  call coop_dictionary_lookup(settings, 'map_maxval', map_maxval, 500.)
  call coop_dictionary_lookup(settings, 'map_minval', map_minval, -500.)


  nmaps = len_trim(map_genre)
  numcls = nmaps*(nmaps+1)/2
  
  allocate(spin(nmaps))
  select case(trim(map_genre))
  case("I")
    spin = 0
    fieldnames = "T"
 case("QU")
    spin = (/ 2, 2 /)
    fieldnames = "EB"
 case("IQU")
    spin = (/ 0, 2, 2 /)
    fieldnames = "TEB"
 case default
    write(*,*) "map_genre = "//trim(map_genre)
    stop "XFASTER only supports map_genre = I/QU/IQU"
 end select
     
  call coop_dictionary_lookup(settings, 'map_unit', map_unit, 'muK')

  call coop_dictionary_lookup(settings, 'num_channels', num_channels)
  numcross = num_channels*(num_channels+1)/2
  matdim = num_channels*nmaps
  allocate(channel_name(num_channels))
  do ich = 1, num_channels
     call coop_dictionary_lookup(settings, 'channel'//COOP_STR_OF(ich)//'_name', channel_name(ich), COOP_STR_OF(ich))
     if(feedback .ge. 2)write(*,*) "channel : "//CHIND(ich)
  enddo
  
  call coop_dictionary_lookup(settings, 'sim_index_width', sim_index_width, 0)

  call coop_dictionary_lookup(settings, "num_ell_bins", num_ell_bins)
  num_qbs = numcross*numcls*num_ell_bins
  allocate(qb(num_qbs))

  num_mubs = numcls*num_ell_bins
  allocate(mub(num_mubs))

  qb = 1.d0
  mub = 1.d0
  call coop_dictionary_lookup(settings, 'lmin_data', lmin_data, 2)
  call coop_dictionary_lookup(settings, 'lmax_data', lmax_data)
  coop_healpix_want_cls = .true.
  if(num_channels .lt. 1)then
     write(*,*) "number of channels = ", num_channels
     stop "XFASTER only supports num_channels >=1"
  endif
  select case(trim(action))
  case("DO_MASK")
     call do_mask()
  case("DO_KERNEL")
     call do_kernel()
  case("DO_NOISE")
     call do_noise()
  case("DO_SIGNAL")
     call do_signal()  
  case("DO_DATA")
     call do_data()
  case("DO_QB")
     do iter = 1, num_iterations
        call do_qb(iter)  !!compute q_b's
     enddo
  case("DO_MUB")
     do iter = 1, num_iterations
        call do_mub(iter, want_plot = (iter .eq. num_iterations))  !!compute mu_b's
     enddo
  case("DO_ALL")
     call do_mask()
     call do_kernel()
     call do_noise()
     call do_signal()
     call do_data()
     do iter = 1, num_iterations
        call do_qb(iter)
     enddo
     do iter = 1, num_iterations
        call do_mub(iter, want_plot = (iter .eq. num_iterations))
     enddo
  case("DO_ALL_BUT_KERNEL")
     call do_noise()
     call do_signal()
     call do_data()
     do iter = 1, num_iterations
        call do_qb(iter)
     enddo
     do iter = 1, num_iterations
        call do_mub(iter, want_plot = (iter .eq. num_iterations))
     enddo
  case("UPDATE_DATA")
     call do_data()
     do iter = 1, num_iterations
        call do_qb(iter)
     enddo
     do iter = 1, num_iterations
        call do_mub(iter, want_plot = (iter .eq. num_iterations))
     enddo
  case("COMPARE_PSEUDO_CLS", "COMPARE_CLS")
     call compare_pseudo_cls()
  case default
     write(*,*) "action  = "//trim(action)
     stop "Unknown action"
  end select

contains

  subroutine do_mask()
    COOP_STRING:: cond_root, hits_root, mask_root, data_map_root
    type(coop_healpix_maps)::mask, cond, hits, premask
    COOP_SINGLE::condmin, hitsmax, maskcut
    COOP_REAL::RA_min, RA_max, DEC_min, DEC_max, theta_min, theta_max, theta, phi, smooth_mask_fwhm
    COOP_INT::ich, nside, pix
    logical::do_RADEC_cut, has_cond, has_hits
    nside = 0
    call coop_dictionary_lookup(settings, "mask_root", mask_root)    
    call coop_dictionary_lookup(settings, "hits_root", hits_root, "")
    if(trim(hits_root).ne. "")then
       call hits%read(trim(hits_root)//"_"//CHIND(1)//".fits", nmaps_wanted = 0)
       nside = hits%nside
       call coop_dictionary_lookup(settings, "hits_cut", hitsmax, -1.e30)
       has_hits = (hitsmax .gt. 0.)
    else
       has_hits = .false.
    endif

    call coop_dictionary_lookup(settings, "cond_root", cond_root,"")
    if(trim(cond_root).ne."")then
       if(nside .eq. 0)then
          call cond%read(trim(cond_root)//"_"//CHIND(1)//".fits", nmaps_wanted = 0)
          nside = cond%nside
       endif
       call coop_dictionary_lookup(settings, "cond_cut", condmin, -1.e30)
       has_cond  = (condmin .gt. 0.)
    else
       has_cond = .false.
    endif

    if(nside .eq. 0)then
       call coop_dictionary_lookup(settings, "data_map_root", data_map_root, "")
       if(trim(data_map_root) .ne. "")then
          call hits%read(trim(data_map_root)//"_"//CHIND(1)//".fits", nmaps_wanted = 0)
          nside = hits%nside
       else
          call coop_dictionary_lookup(settings, "nside", nside)
       endif
    endif
    call coop_dictionary_lookup(settings, "RA_min", RA_min, 1.d30)
    call coop_dictionary_lookup(settings, "RA_max", RA_max, -1.d30)
    call coop_dictionary_lookup(settings, "DEC_min", DEC_min, 1.d30)
    call coop_dictionary_lookup(settings, "DEC_max", DEC_max, -1.d30)
    do_RADEC_cut = (abs(RA_min) .lt. 1.d20 .and. abs(RA_max) .lt. 1.d20 .and. abs(dec_min) .lt. 1.d20 .and. abs(dec_max) .lt. 1.d20)
    if(do_RADEC_cut)then
       ra_min = ra_min * coop_SI_degree
       ra_max = ra_max * coop_SI_degree
       dec_min = dec_min * coop_SI_degree
       dec_max = dec_max * coop_SI_degree
       call coop_periodic_select_range(ra_min, -coop_pi, coop_pi)
       call coop_periodic_select_range(ra_max, ra_min, ra_min + coop_2pi)
       call coop_periodic_select_range(dec_min, -coop_pio2, coop_pio2)
       call coop_periodic_select_range(dec_max, -coop_pio2, coop_pio2)
       theta_min = coop_pio2 - max(dec_min, dec_max)
       theta_max = coop_pio2 - min(dec_min, dec_max)
       call premask%init(nmaps = 1, nside = nside, genre="MASK")
       !$omp parallel do private(pix, theta, phi)
       do pix = 0, premask%npix-1
          call premask%pix2ang(pix, theta, phi)
          call coop_periodic_select_range(phi, ra_min, ra_min + coop_2pi)
          if(theta .lt. theta_min .or. theta.gt. theta_max .or. phi .gt. ra_max)then
             premask%map(pix, 1) = 0.
          else
             premask%map(pix, 1) = 1.
          endif
       enddo
       !$omp end parallel do

    endif
    call coop_dictionary_lookup(settings, "smooth_mask_arcmin", smooth_mask_fwhm, 0.d0)
    smooth_mask_fwhm = smooth_mask_fwhm * coop_SI_arcmin
    do ich = 1, num_channels
       if(has_cond)call cond%read(trim(cond_root)//"_"//CHIND(ich)//".fits", nmaps_wanted = 1)
       if(has_hits)call hits%read(trim(hits_root)//"_"//CHIND(ich)//".fits", nmaps_wanted = 1)
       if(do_radec_cut)then
          mask = premask
       else
          call mask%init(nmaps = 1, nside = nside, genre="MASK")
          mask%map(:,1) = 1.
       endif
       if(has_hits)then
          where(hits%map(:,1) .lt. hitsmax)
             mask%map(:,1) = 0.
          end where
       end if
       if(has_cond)then
          where(cond%map(:,1) .gt. condmin)
             mask%map(:,1) = 0.
          end where
       endif
       if(smooth_mask_fwhm .gt. 0.d0)then
          call mask%smooth(fwhm = smooth_mask_fwhm)
       endif
       if(feedback .ge. 1) write(*,*) "fsky for mask #"//COOP_STR_OF(ich)//" is "//COOP_STR_OF(sum(mask%map(:,1))/mask%npix)
       call mask%write(trim(mask_root)//"_"//CHIND(ich)//".fits")
       if(feedback .ge. 1) write(*,*) "channel "//COOP_STR_OF(ich)//" fsky: "//COOP_STR_OF(sum(mask%map(:,1))/mask%npix)
    enddo
  end subroutine do_mask

  subroutine do_kernel()
    COOP_STRING::kernel_root, mask_root
    COOP_INT::ich1, ich2, lmax_mask, lmax_kernel
    type(coop_healpix_maps)::mask1, mask2
    type(coop_dictionary)::header
    COOP_REAL,dimension(:,:,:),allocatable::kernel
    if(feedback .ge. 1)then
       write(*,*) "***************************************************"
       write(*,*) "Doing mask kernels."
    endif
    call coop_dictionary_lookup(settings, "kernel_root", kernel_root)
    call coop_dictionary_lookup(settings, "mask_root", mask_root)
    call coop_dictionary_lookup(settings, "lmax_mask", lmax_mask, 300)
    call coop_dictionary_lookup(settings, "lmax_kernel", lmax_kernel, lmax_data)
    allocate(kernel(0:lmax_kernel, 0:lmax_kernel, 4))
    call header%insert("LMIN", "0")
    call header%insert("LMAX", COOP_STR_OF(lmax_kernel))
    do ich1=1, num_channels
       call mask1%read(trim(mask_root)//"_"//CHIND(ich1)//".fits", nmaps_wanted=1, nested = .false.)
       call mask1%map2alm(lmax = lmax_mask)
       call header%insert("WEIGHT", COOP_STR_OF( sum(mask1%map(:,1)**2)**2/dble(mask1%npix)/sum(mask1%map(:,1)**4)), overwrite = .true.)
       call coop_pseudoCl_get_kernel(lmax_mask = lmax_mask, Cl_mask = dble(mask1%cl(0:lmax_mask, 1)), lmin = 0 , lmax= lmax_kernel, kernel = kernel)
       call coop_fits_file_write_image_3d(filename = trim(kernel_root)//"_"//CHIND(ich1)//"_"//CHIND(ich1)//".fits", image = kernel, header=header)
       if(feedback .ge. 1) &
            write(*,*) "kernel file "//trim(kernel_root)//"_"//CHIND(ich1)//"_"//CHIND(ich1)//".fits is produced."
       do ich2 = ich1+1, num_channels
          call mask2%read(trim(mask_root)//"_"//CHIND(ich2)//".fits", nmaps_wanted=1, nested = .false.)
          call mask2%map2alm(lmax = lmax_mask)
          call mask2%get_cls(mask1)
          call coop_pseudoCl_get_kernel(lmax_mask = lmax_mask, Cl_mask = dble(mask2%cl(0:lmax_mask, 1)), lmin = 0 , lmax= lmax_kernel, kernel = kernel)
          call header%insert("WEIGHT", COOP_STR_OF(sum(mask1%map(:,1)*mask2%map(:,1))**2/mask1%npix/sum((mask1%map(:,1)*mask2%map(:,1))**2)), overwrite = .true.)
          call coop_fits_file_write_image_3d(filename = trim(kernel_root)//"_"//CHIND(ich1)//"_"//CHIND(ich2)//".fits", image = kernel, header = header)
          if(feedback .ge. 1) &
               write(*,*) "kernel file "//trim(kernel_root)//"_"//CHIND(ich1)//"_"//CHIND(ich2)//".fits is produced."
       enddo
    enddo
    call mask1%free()
    call mask2%free()
    deallocate(kernel)
  end subroutine do_kernel

  subroutine do_noise()
    COOP_STRING::noise_sim_root, mask_root, noise_cl_root
    type(coop_healpix_maps)::noise1, noise2
    type(coop_cls)::noisepower
    COOP_INT::ich1, ich2, isim, num_noise_sims
    call coop_dictionary_lookup(settings, "noise_sim_root", noise_sim_root)
    call coop_dictionary_lookup(settings, "noise_cl_root", noise_cl_root)
    call coop_dictionary_lookup(settings, "mask_root", mask_root)
    call coop_dictionary_lookup(settings, "num_noise_sims", num_noise_sims)
    call noisepower%init(lmin = lmin_data, lmax = lmax_data, genre = map_genre, unit = map_unit, spin = spin)
    if(feedback .ge. 1)then
       write(*,*) "***************************************************"
       write(*,*) "computing noise power from "//COOP_STR_OF(num_noise_sims)//" simulations"
    endif
    do ich1 = 1, num_channels
       do ich2 = ich1, num_channels
          if(feedback .ge. 1) &
               write(*,"(A$)") "channel "//CHIND(ich1)//" x channel "//CHIND(ich2)
          noisepower%cls = 0.d0
          do isim = 1, num_noise_sims
             if(feedback .ge. 1) &
                  write(*,"(A$)") "."
             call read_map(noise1, trim(noise_sim_root)//"_"//CHSIMIND(ich1, isim)//".fits", trim(mask_root)//"_"//CHIND(ich1)//".fits")
             if(ich2.ne.ich1)then
                call read_map(noise2, trim(noise_sim_root)//"_"//CHSIMIND(ich2, isim)//".fits", trim(mask_root)//"_"//CHIND(ich2)//".fits")
                call noise1%get_cls(noise2)
             endif
             noisepower%cls = noisepower%cls + noise1%cl(lmin_data:lmax_data, :)
          enddo
          noisepower%cls = noisepower%cls/num_noise_sims
          call noisepower%smooth(5)
          call noisepower%dump(trim(noise_cl_root)//"_"//CHIND(ich1)//"_"//CHIND(ich2)//".fits")
          if(feedback .ge. 1) &
               write(*,*)
       enddo
    enddo
    call noise1%free()
    call noise2%free()
    call noisepower%free()
  end subroutine do_noise


  subroutine do_signal()
    COOP_STRING::signal_sim_root, signal_cl_root, mask_root
    logical want_T, want_EB
    type(coop_cls)::signalpower
    COOP_INT::ich1, ich2, num_signal_sims, isim, ib
    type(coop_healpix_maps)::sig1, sig2
    COOP_REAL::kernel(lmin_data:lmax_data, lmin_data:lmax_data, 4), weight
    COOP_REAL::Cl_pseudo(lmin_data:lmax_data, 6), Cl(lmin_data:lmax_data, 6)
    call coop_dictionary_lookup(settings, "mask_root", mask_root)
    call coop_dictionary_lookup(settings, "signal_sim_root", signal_sim_root)
    call coop_dictionary_lookup(settings, "signal_cl_root", signal_cl_root)
    call coop_dictionary_lookup(settings, "num_signal_sims", num_signal_sims)
    if(feedback .ge. 1)then
       write(*,*) "***************************************************"
       write(*,*) "computing sigal power from "//COOP_STR_OF(num_signal_sims)//" simulations"
    endif
    call signalpower%init(lmin = lmin_data, lmax = lmax_data, genre = map_genre, unit = map_unit, spin = spin)
    do ich1 = 1, num_channels
       do ich2 = ich1, num_channels
          if(feedback .ge. 1) &
               write(*,"(A$)") "channel "//CHIND(ich1)//" x channel "//CHIND(ich2)
          signalpower%cls = 0.d0
          do isim = 1, num_signal_sims
             if(feedback .ge. 1) &
                  write(*,"(A$)") "."
             call read_map(sig1, trim(signal_sim_root)//"_"//CHSIMIND(ich1, isim)//".fits", trim(mask_root)//"_"//CHIND(ich1)//".fits")
             if(ich2.ne.ich1)then
                call read_map(sig2, trim(signal_sim_root)//"_"//CHSIMIND(ich2, isim)//".fits", trim(mask_root)//"_"//CHIND(ich2)//".fits")
                call sig1%get_cls(sig2)
             endif
             signalpower%cls = signalpower%cls + sig1%cl(lmin_data:lmax_data, :)
          enddo
          signalpower%cls = signalpower%cls/num_signal_sims
          call signalpower%dump(trim(signal_cl_root)//"_"//CHIND(ich1)//"_"//CHIND(ich2)//".fits")
          if(feedback .ge. 1) &
               write(*,*)
       enddo
    enddo
    call sig1%free()
    call sig2%free()
    call signalpower%free()
  end subroutine do_signal


  subroutine do_data()
    COOP_STRING::data_map_root, data_cl_root, mask_root
    type(coop_healpix_maps)::map1, map2
    COOP_INT::ich1, ich2
    type(coop_cls)::datapower
    if(feedback .ge. 1) then
       write(*,*) "***************************************************"
       write(*,"(A$)") "Doing pseudo Cls"
    endif
    call coop_dictionary_lookup(settings, "data_map_root", data_map_root)
    call coop_dictionary_lookup(settings, "data_cl_root", data_cl_root)
    call coop_dictionary_lookup(settings, "mask_root", mask_root)
    do ich1 =1, num_channels
       call read_map(map1, trim(data_map_root)//"_"//CHIND(ich1)//".fits", trim(mask_root)//"_"//CHIND(ich1)//".fits")
       call datapower%init(lmin = lmin_data, lmax = lmax_data, genre = map_genre, unit = map_unit, spin = spin, cls = dble(map1%cl(lmin_data:lmax_data, :)))
       call datapower%dump(trim(data_cl_root)//"_"//CHIND(ich1)//"_"//CHIND(ich1)//".fits")
       if(feedback .ge.1 ) &
            write(*,"(A$)") "."
       do ich2 = ich1+1, num_channels
          call read_map(map2, trim(data_map_root)//"_"//CHIND(ich2)//".fits", trim(mask_root)//"_"//CHIND(ich2)//".fits")
          call map2%get_cls(map1)
          call datapower%init(lmin = lmin_data, lmax = lmax_data, genre = map_genre, unit = map_unit, spin = spin, cls = dble(map2%cl(lmin_data:lmax_data, :)))
          call datapower%dump(trim(data_cl_root)//"_"//CHIND(ich1)//"_"//CHIND(ich2)//".fits")
          if(feedback .ge.1 ) &
               write(*,"(A$)") "."
       enddo
    enddo
    if(feedback .ge. 1) &
         write(*,*)
    call map1%free()
    call map2%free()
    call datapower%free()
  end subroutine do_data
  

  subroutine do_qb(iter)
    COOP_STRING:: kernel_root, signal_cl_root, noise_cl_root
    COOP_INT:: l, iter, num_qb_used, info
    type(coop_cls),dimension(numcross)::signalpower, noisepower
    type(coop_binned_cls)::templatepower
    COOP_REAL::kernel(lmin_data:lmax_data, lmin_data:lmax_data, 4,numcross)
    COOP_REAL::weight(numcross)
    COOP_INT::i, j, ind, ii, jj, ib, ich1, ich2, icl, jb, jch1, jch2, jcl, iqb, ibase1, ibase2
    COOP_REAL::mat(matdim, matdim), dSdq(matdim, matdim, num_qbs), wmat(matdim, matdim), datamat(matdim, matdim), invDdSdq(matdim, matdim, num_qbs)
    logical::nonzero(num_qbs)
    COOP_REAL::Fisher(num_qbs, num_qbs)
    COOP_REAL::vecb(num_qbs)
    COOP_REAL::Cb_EE(numcls, num_ell_bins), Cb_BB(numcls, num_ell_bins)
    type(coop_file)::fp
    COOP_REAL,dimension(:,:),allocatable::Fisher_used
    COOP_INT,dimension(:),allocatable::index_qb_used
    if(feedback .ge. 1)then
       write(*,*) "***************************************************"
       write(*,*) "Iterating the maximum-likelihood Cls"
    endif
    call coop_dictionary_lookup(settings, "kernel_root", kernel_root)
    call coop_dictionary_lookup(settings, "signal_cl_root", signal_cl_root)
    call coop_dictionary_lookup(settings, "noise_cl_root", noise_cl_root)
    call templatepower%load(model_cl_file)
    call templatepower%filter(lmin = lmin_data, lmax = lmax_data)
    call templatepower%alloc(num_ell_bins)

    do i = 1, num_channels
       do j = i, num_channels
          ind = COOP_MATSYM_INDEX(num_channels, i,j)
          call load_kernel(trim(kernel_root)//"_"//CHIND(i)//"_"//CHIND(j)//".fits", lmin_data, lmax_data, kernel(:,:,:,ind), weight(ind))
          call signalpower(ind)%load(trim(signal_cl_root)//"_"//CHIND(i)//"_"//CHIND(j)//".fits")
          call signalpower(ind)%filter(lmin = lmin_data, lmax = lmax_data)
          call noisepower(ind)%load(trim(noise_cl_root)//"_"//CHIND(i)//"_"//CHIND(j)//".fits")
          call noisepower(ind)%filter(lmin = lmin_data, lmax = lmax_data)
       enddo
    enddo
    vecb = 0.d0
    Fisher = 0.d0
    do l = lmin_data, lmax_data
       dSdq = 0.d0
       do ich1 = 1, num_channels
          ibase1 = (ich1-1) * nmaps
          do ich2 = ich1, num_channels
             ind = COOP_MATSYM_INDEX(num_channels, ich1, ich2)
             ibase2 = (ich2-1) * nmaps
             call get_Cbs(templatepower, l, kernel(:,:,:,ind), Cb_EE, Cb_BB)
             select case(trim(map_genre))
             case("I")
                do ib = 1, num_ell_bins
                   dSdq(ich1, ich2, QB_IND(ib, ind, 1)) = Cb_EE(1, ib)
                enddo
             case("QU")
                do ib = 1, num_ell_bins
                   dSdq(ibase1+coop_EB_index_E, ibase2+coop_EB_index_E, QB_IND(ib, ind, coop_EB_index_EE)) = Cb_EE(coop_EB_index_EE, ib)
                   dSdq(ibase1+coop_EB_index_B, ibase2+coop_EB_index_B, QB_IND(ib, ind, coop_EB_index_EE)) = Cb_EE(coop_EB_index_BB, ib)

                   dSdq(ibase1+coop_EB_index_E, ibase2+coop_EB_index_E, QB_IND(ib, ind, coop_EB_index_BB)) = Cb_BB(coop_EB_index_EE, ib)
                   dSdq(ibase1+coop_EB_index_B, ibase2+coop_EB_index_B, QB_IND(ib, ind, coop_EB_index_BB)) = Cb_BB(coop_EB_index_BB, ib)

                   dSdq(ibase1+coop_EB_index_E, ibase2+coop_EB_index_B, QB_IND(ib, ind, coop_EB_index_EB)) = Cb_EE(coop_EB_index_EB, ib)
                   dSdq(ibase1+coop_EB_index_B, ibase2+coop_EB_index_E, QB_IND(ib, ind, coop_EB_index_EB)) = Cb_EE(coop_EB_index_EB, ib)
                   
                enddo
             case("IQU")
                do ib = 1, num_ell_bins
                   dSdq(ibase1+coop_TEB_index_E, ibase2+coop_TEB_index_E, QB_IND(ib, ind, coop_TEB_index_EE)) = Cb_EE(coop_TEB_index_EE, ib)
                   dSdq(ibase1+coop_TEB_index_B, ibase2+coop_TEB_index_B, QB_IND(ib, ind, coop_TEB_index_EE)) = Cb_EE(coop_TEB_index_BB, ib)

                   dSdq(ibase1+coop_TEB_index_E, ibase2+coop_TEB_index_E, QB_IND(ib, ind, coop_TEB_index_BB)) = Cb_BB(coop_TEB_index_EE, ib)
                   dSdq(ibase1+coop_TEB_index_B, ibase2+coop_TEB_index_B, QB_IND(ib, ind, coop_TEB_index_BB)) = Cb_BB(coop_TEB_index_BB, ib)

                   dSdq(ibase1+coop_TEB_index_E, ibase2+coop_TEB_index_B, QB_IND(ib, ind, coop_TEB_index_EB)) = Cb_EE(coop_TEB_index_EB, ib)
                   dSdq(ibase1+coop_TEB_index_B, ibase2+coop_TEB_index_E, QB_IND(ib, ind, coop_TEB_index_EB)) = Cb_EE(coop_TEB_index_EB, ib)

                   dSdq(ibase1+coop_TEB_index_T, ibase2+coop_TEB_index_T, QB_IND(ib, ind, coop_TEB_index_TT)) = Cb_EE(coop_TEB_index_TT, ib)
                   dSdq(ibase1+coop_TEB_index_T, ibase2+coop_TEB_index_E, QB_IND(ib, ind, coop_TEB_index_TE)) = Cb_EE(coop_TEB_index_TE, ib)
                   dSdq(ibase1+coop_TEB_index_T, ibase2+coop_TEB_index_B, QB_IND(ib, ind, coop_TEB_index_TB)) = Cb_EE(coop_TEB_index_TB, ib)

                   dSdq(ibase1+coop_TEB_index_E, ibase2+coop_TEB_index_T, QB_IND(ib, ind, coop_TEB_index_TE)) = Cb_EE(coop_TEB_index_TE, ib)
                   dSdq(ibase1+coop_TEB_index_B, ibase2+coop_TEB_index_T, QB_IND(ib, ind, coop_TEB_index_TB)) = Cb_EE(coop_TEB_index_TB, ib)
                enddo
             case default
                write(*,*) "map_genre = "//trim(map_genre)
                stop "XFASTER only supports map_genre = I/QU/IQU"
             end select

             do i = 1, nmaps
                do j = i, nmaps
                   mat(ibase1 + i, ibase2+j) = noisepower(ind)%cls(l, COOP_MATSYM_INDEX(nmaps, i, j))
                   datamat(ibase1 + i, ibase2+j) = signalpower(ind)%cls(l, COOP_MATSYM_INDEX(nmaps, i, j)) 
                   if(i.ne.j)then
                      mat(ibase1 + j, ibase2+i) = mat(ibase1 + i, ibase2 + j)
                      datamat(ibase1 + j, ibase2+i) = datamat(ibase1 + i, ibase2 + j)
                   endif

                enddo
             enddo
             if(ich1 .ne. ich2)then
                dSdq(ibase2+1:ich2*nmaps, ibase1+1:ich1*nmaps, :) = dSdq(ibase1+1:ich1*nmaps, ibase2+1:ich2*nmaps, :)
                mat(ibase2+1:ich2*nmaps, ibase1+1:ich1*nmaps) = mat(ibase1+1:ich1*nmaps, ibase2+1:ich2*nmaps)
                datamat(ibase2+1:ich2*nmaps, ibase1+1:ich1*nmaps) = datamat(ibase1+1:ich1*nmaps, ibase2+1:ich2*nmaps)
             endif
          enddo
       enddo
       do ib = 1, num_qbs
          if(all( dSdq(:,:,ib) .eq. 0.d0))then
             nonzero(ib) = .false.
          else
             mat = mat + dSdq(:,:, ib)*qb(ib)
             nonzero(ib) = .true.
          endif
       enddo
       call coop_sympos_inverse(matdim, matdim, mat,info)
       if(info .ne. 0)then
          stop "do_qb does not converge, you may try reducing the number of bins."
       endif
       do ich1 = 1, num_channels
          ibase1 = (ich1-1) * nmaps
          do ich2 = 1, num_channels
             ibase2 = (ich2-1) * nmaps
             wmat(ibase1+1:ich1*nmaps, ibase2+1:ich2*nmaps) = mat(ibase1+1:ich1*nmaps, ibase2+1:ich2*nmaps) *weight(COOP_MATSYM_INDEX(num_channels, ich1, ich2))
          enddo
       enddo
       !$omp parallel do
       do ib = 1, num_qbs
          if(nonzero(ib))then
             invDdSdq(:, :, ib) = matmul(wmat, matmul(dSdq(:,:, ib), mat))
          endif
       enddo
       !$omp end parallel do
       !$omp parallel do private(ib, jb)
       do ib = 1, num_qbs
          if(nonzero(ib))then
             do jb = ib, num_qbs
                if(nonzero(jb))then
                   Fisher(ib, jb) = Fisher(ib, jb) + coop_matrix_product_trace(invDdSdq(:,:,ib), dSdq(:,:,jb)) * (l+0.5d0)
                endif
             enddo
             vecb(ib) = vecb(ib) + coop_matrix_product_trace(invDdSdq(:,:,ib),  datamat)*(l+0.5d0)
          endif
       enddo
       !$omp end parallel do
    enddo
    if(feedback .ge. 1)then
       write(*,*)
    endif
    do ib = 1, num_qbs
       do jb = ib+1, num_qbs
          Fisher(jb, ib)  = Fisher(ib, jb)
       enddo
    enddo
    do ib=1, num_qbs
       nonzero(ib) = Fisher(ib, ib) .gt. fisher_threshold
    enddo
    num_qb_used = count(nonzero)
    allocate(index_qb_used(num_qb_used), Fisher_used(num_qb_used, num_qb_used))
    i = 0
    do ib = 1, num_qbs
       if(nonzero(ib))then
          i = i + 1
          index_qb_used(i) = ib
       endif
    enddo
    Fisher_used = Fisher(index_qb_used, index_qb_used)
    call coop_sympos_inverse(num_qb_used, num_qb_used, Fisher_used)
    qb(index_qb_used) = matmul(Fisher_used, vecb(index_qb_used))
    call fp%open(trim(qb_output_root)//"_ITER"//COOP_STR_OF(iter)//".dat")
    do ib = 1, num_ell_bins
       write(fp%unit, "(F10.2, "//COOP_STR_OF(numcross*numcls)//"E16.7)") templatepower%lb(ib), qb((ib-1)*numcross*numcls+1:ib*numcross*numcls)
    enddo
    call fp%close()
    do ind = 1, numcross
       call noisepower(ind)%free()
       call signalpower(ind)%free()
    enddo
    call templatepower%free()
    deallocate(index_qb_used, Fisher_used)
  end subroutine do_qb


  subroutine do_mub(iter, want_plot)
    logical,optional::want_plot
    COOP_STRING:: kernel_root, data_cl_root, noise_cl_root
    COOP_INT:: l, iter, num_mub_used
    type(coop_cls),dimension(numcross)::signalpower, noisepower
    type(coop_binned_cls)::templatepower
    type(coop_cls)::modelcl
    COOP_REAL::kernel(lmin_data:lmax_data, lmin_data:lmax_data, 4,numcross)
    COOP_REAL::weight(numcross)
    COOP_INT::i, j, ind, ii, jj, ib, ich1, ich2, icl, jb, jch1, jch2, jcl, iqb, mapind, ibase1, ibase2, info
    COOP_REAL::mat(matdim, matdim), dSdq(matdim, matdim, num_mubs), wmat(matdim, matdim), datamat(matdim, matdim), invDdSdq(matdim, matdim, num_mubs), ells(lmin_data:lmax_data), tr
    logical::nonzero(num_mubs)
    COOP_REAL::Fisher(num_mubs, num_mubs)
    COOP_REAL::vecb(num_mubs), errorbars(num_mubs), lambda
    COOP_REAL::Cb_EE(numcls, num_ell_bins), Cb_BB(numcls, num_ell_bins)
    type(coop_file)::fp
    COOP_REAL,dimension(:,:),allocatable::Fisher_used
    COOP_INT,dimension(:),allocatable::index_mub_used
    type(coop_asy)::figure
    if(feedback .ge. 1)then
       write(*,*) "***************************************************"
       write(*,*) "Iterating the maximum-likelihood Cls"
    endif
    call coop_dictionary_lookup(settings, "kernel_root", kernel_root)
    call coop_dictionary_lookup(settings, "data_cl_root", data_cl_root)
    call coop_dictionary_lookup(settings, "noise_cl_root", noise_cl_root)
    call load_qb(num_iterations)  !!get the converged q_b's

    call templatepower%load(model_cl_file)
    call templatepower%filter(lmin = lmin_data, lmax = lmax_data)
    call templatepower%alloc(num_ell_bins)

    do i = 1, num_channels
       do j = i, num_channels
          ind = COOP_MATSYM_INDEX(num_channels, i,j)
          call load_kernel(trim(kernel_root)//"_"//CHIND(i)//"_"//CHIND(j)//".fits", lmin_data, lmax_data, kernel(:,:,:,ind), weight(ind))
          call signalpower(ind)%load(trim(data_cl_root)//"_"//CHIND(i)//"_"//CHIND(j)//".fits")
          call signalpower(ind)%filter(lmin = lmin_data, lmax = lmax_data)
          call noisepower(ind)%load(trim(noise_cl_root)//"_"//CHIND(i)//"_"//CHIND(j)//".fits")
          call noisepower(ind)%filter(lmin = lmin_data, lmax = lmax_data)
          signalpower(ind)%cls = signalpower(ind)%cls - noisepower(ind)%cls !!remove noise
       enddo
    enddo
    vecb = 0.d0
    Fisher = 0.d0
    do l = lmin_data, lmax_data
       dSdq = 0.d0
       do ich1 = 1, num_channels
          ibase1 = (ich1-1) * nmaps
          do ich2 = ich1, num_channels
             ind = COOP_MATSYM_INDEX(num_channels, ich1, ich2)
             ibase2 = (ich2-1) * nmaps
             call get_Cbs(templatepower, l, kernel(:,:,:,ind), Cb_EE, Cb_BB)
             select case(trim(map_genre))
             case("I")
                do ib = 1, num_ell_bins
                   dSdq(ich1, ich2, MUB_IND(ib, 1)) =   Cb_EE(1, ib)*qb(QB_IND(ib, ind, 1))
                enddo
             case("QU")
                do ib = 1, num_ell_bins
                   dSdq(ibase1+coop_EB_index_E, ibase2+coop_EB_index_E, MUB_IND(ib, coop_EB_index_EE)) =  Cb_EE(coop_EB_index_EE, ib) * qb(QB_IND(ib, ind, coop_EB_index_EE))
                   dSdq(ibase1+coop_EB_index_B, ibase2+coop_EB_index_B, MUB_IND(ib, coop_EB_index_EE)) =  Cb_EE(coop_EB_index_BB, ib) * qb(QB_IND(ib, ind, coop_EB_index_EE))
                   dSdq(ibase1+coop_EB_index_E, ibase2+coop_EB_index_E, MUB_IND(ib, coop_EB_index_BB)) =  Cb_BB(coop_EB_index_EE, ib) * qb(QB_IND(ib, ind, coop_EB_index_BB))
                   dSdq(ibase1+coop_EB_index_B, ibase2+coop_EB_index_B, MUB_IND(ib, coop_EB_index_BB)) =  Cb_BB(coop_EB_index_BB, ib) * qb(QB_IND(ib, ind, coop_EB_index_BB))
                   dSdq(ibase1+coop_EB_index_E, ibase2+coop_EB_index_B, MUB_IND(ib, coop_EB_index_EB)) =  Cb_EE(coop_EB_index_EB, ib) * qb(QB_IND(ib, ind, coop_EB_index_EB))
                   dSdq(ibase1+coop_EB_index_B, ibase2+coop_EB_index_E, MUB_IND(ib, coop_EB_index_EB)) =  dSdq(ibase1+coop_EB_index_E, ibase2+coop_EB_index_B, MUB_IND(ib, coop_EB_index_EB))
                   
                enddo
             case("IQU")
                do ib = 1, num_ell_bins
                   dSdq(ibase1+coop_TEB_index_E, ibase2+coop_TEB_index_E, MUB_IND(ib, coop_TEB_index_EE)) = Cb_EE(coop_TEB_index_EE, ib)*qb(QB_IND(ib, ind, coop_TEB_index_EE))
                   dSdq(ibase1+coop_TEB_index_B, ibase2+coop_TEB_index_B, MUB_IND(ib, coop_TEB_index_EE)) = Cb_EE(coop_TEB_index_BB, ib)*qb(QB_IND(ib, ind, coop_TEB_index_EE))

                   dSdq(ibase1+coop_TEB_index_E, ibase2+coop_TEB_index_E, MUB_IND(ib,coop_TEB_index_BB)) =  Cb_BB(coop_TEB_index_EE, ib)*qb(QB_IND(ib, ind, coop_TEB_index_BB))
                   dSdq(ibase1+coop_TEB_index_B, ibase2+coop_TEB_index_B, MUB_IND(ib, coop_TEB_index_BB)) = Cb_BB(coop_TEB_index_BB, ib)*qb(QB_IND(ib, ind, coop_TEB_index_BB))

                   dSdq(ibase1+coop_TEB_index_E, ibase2+coop_TEB_index_B, MUB_IND(ib, coop_TEB_index_EB)) = Cb_EE(coop_TEB_index_EB, ib)*qb(QB_IND(ib, ind, coop_TEB_index_EB))
                   dSdq(ibase1+coop_TEB_index_B, ibase2+coop_TEB_index_E, MUB_IND(ib, coop_TEB_index_EB)) =  dSdq(ibase1+coop_TEB_index_E, ibase2+coop_TEB_index_B, MUB_IND(ib, coop_TEB_index_EB))

                   dSdq(ibase1+coop_TEB_index_T, ibase2+coop_TEB_index_T, MUB_IND(ib, coop_TEB_index_TT)) = Cb_EE(coop_TEB_index_TT, ib)*qb(QB_IND(ib, ind, coop_TEB_index_TT))

                   dSdq(ibase1+coop_TEB_index_T, ibase2+coop_TEB_index_E, MUB_IND(ib, coop_TEB_index_TE)) = Cb_EE(coop_TEB_index_TE, ib)*qb(QB_IND(ib, ind, coop_TEB_index_TE))
                   dSdq(ibase1+coop_TEB_index_T, ibase2+coop_TEB_index_B, MUB_IND(ib, coop_TEB_index_TB)) = Cb_EE(coop_TEB_index_TB, ib)*qb(QB_IND(ib, ind, coop_TEB_index_TB))

                   dSdq(ibase1+coop_TEB_index_E, ibase2+coop_TEB_index_T, MUB_IND(ib, coop_TEB_index_TE)) =  dSdq(ibase1+coop_TEB_index_T, ibase2+coop_TEB_index_E, MUB_IND(ib, coop_TEB_index_TE))
                   dSdq(ibase1+coop_TEB_index_B, ibase2+coop_TEB_index_T, MUB_IND(ib, coop_TEB_index_TB)) =  dSdq(ibase1+coop_TEB_index_T, ibase2+coop_TEB_index_B, MUB_IND(ib, coop_TEB_index_TB))
                enddo
             case default
                write(*,*) "map_genre = "//trim(map_genre)
                stop "XFASTER only supports map_genre = I/QU/IQU"
             end select

             do i = 1, nmaps
                do j = i, nmaps
                   mapind = COOP_MATSYM_INDEX(nmaps, i, j)
                   mat(ibase1 + i, ibase2+j) = noisepower(ind)%cls(l, mapind)
                   datamat(ibase1 + i, ibase2+j) = signalpower(ind)%cls(l, mapind) 
                   if(i.ne.j)then
                      mat(ibase1 + j, ibase2+i) = mat(ibase1 + i, ibase2 + j)
                      datamat(ibase1 + j, ibase2+i) = datamat(ibase1 + i, ibase2 + j)
                   endif

                enddo
             enddo
             if(ich1 .ne. ich2)then
                dSdq(ibase2+1:ibase2+nmaps, ibase1+1:ibase1+nmaps, :) = dSdq(ibase1+1:ibase1+nmaps, ibase2+1:ibase2+nmaps, :)
                mat(ibase2+1:ibase2+nmaps, ibase1+1:ibase1+nmaps) = mat(ibase1+1:ibase1+nmaps, ibase2+1:ibase2+nmaps)
                datamat(ibase2+1:ibase2+nmaps, ibase1+1:ibase1+nmaps) = datamat(ibase1+1:ibase1+nmaps, ibase2+1:ibase2+nmaps)
             endif
          enddo
       enddo
       do ib = 1, num_mubs
          if(all( dSdq(:,:,ib) .eq. 0.d0))then
             nonzero(ib) = .false.
          else
             mat = mat + dSdq(:,:, ib)*mub(ib)
             nonzero(ib) = .true.
          endif
       enddo
       wmat = mat
       call coop_sympos_inverse(matdim, matdim, mat, info)
       if(info .ne. 0)then
          stop "do_mub does not converge, you may try reducing the number of bins or improve the quality of the mask."
       endif

       do ich1 = 1, num_channels
          ibase1 = (ich1-1) * nmaps
          do ich2 = 1, num_channels
             ibase2 = (ich2-1) * nmaps
             wmat(ibase1+1:ibase1+nmaps, ibase2+1:ibase2+nmaps) = mat(ibase1+1:ibase1+nmaps, ibase2+1:ibase2+nmaps) *weight(COOP_MATSYM_INDEX(num_channels, ich1, ich2))
          enddo
       enddo
       !$omp parallel do
       do ib = 1, num_mubs
          if(nonzero(ib))then
             invDdSdq(:, :, ib) = matmul(wmat, matmul(dSdq(:,:, ib), mat))
          endif
       enddo
       !$omp end parallel do
       !$omp parallel do private(ib, jb)
       do ib = 1, num_mubs
          if(nonzero(ib))then
             do jb = ib, num_mubs
                if(nonzero(jb))then
                   Fisher(ib, jb) = Fisher(ib, jb) + coop_matrix_product_trace(invDdSdq(:,:,ib), dSdq(:,:,jb)) * (l+0.5d0)
                endif
             enddo
             vecb(ib) = vecb(ib) + coop_matrix_product_trace(invDdSdq(:,:,ib),  datamat)*(l+0.5d0)
          endif
       enddo
       !$omp end parallel do
    enddo
    write(*,*)
    do ib = 1, num_mubs
       do jb = ib+1, num_mubs
          Fisher(jb, ib)  = Fisher(ib, jb)
       enddo
    enddo
    do ib=1, num_mubs
       nonzero(ib) = Fisher(ib, ib) .gt. fisher_threshold
    enddo
    num_mub_used = count(nonzero)
    allocate(index_mub_used(num_mub_used), Fisher_used(num_mub_used, num_mub_used))
    i = 0
    do ib = 1, num_mubs
       if(nonzero(ib))then
          i = i + 1
          index_mub_used(i) = ib
       endif
    enddo
    Fisher_used = Fisher(index_mub_used, index_mub_used)
    call coop_sympos_inverse(num_mub_used, num_mub_used, Fisher_used)
    mub(index_mub_used) = matmul(Fisher_used, vecb(index_mub_used))
    mub(index_mub_used) = max(0.2d0, min(mub(index_mub_used), 5.d0)) !!set limits 
    if(feedback .ge. 2) write(*, "("//COOP_STR_OF(num_mub_used)//"F10.3)") mub(index_mub_used)
    errorbars = 0.d0
    do i = 1, num_mub_used
       errorbars(index_mub_used(i)) = sqrt(Fisher_used(i,i))
    enddo
    call fp%open(trim(mub_output_root)//"_ITER"//COOP_STR_OF(iter)//".dat")
    do ib = 1, num_ell_bins
       write(fp%unit, "(F10.2, "//COOP_STR_OF(numcls)//"E16.7)") templatepower%lb(ib), mub((ib-1)*numcls+1:ib*numcls)
    enddo
    call fp%close()

    call fp%open(trim(cl_output_root)//"_BINNED_ITER"//COOP_STR_OF(iter)//".dat")
    do ib = 1, num_ell_bins
       write(fp%unit, "(F10.2, "//COOP_STR_OF(numcls*2)//"E16.7)") templatepower%lb(ib), mub((ib-1)*numcls+1:ib*numcls)*templatepower%Cls(nint(templatepower%lb(ib)),1:numcls)*(templatepower%lb(ib)*( templatepower%lb(ib) + 1.d0)/coop_2pi), abs(errorbars((ib-1)*numcls+1:ib*numcls)*templatepower%Cls(nint(templatepower%lb(ib)),1:numcls))*(templatepower%lb(ib)*( templatepower%lb(ib) + 1.d0)/coop_2pi)
    enddo
    call fp%close()
    if(present(want_plot))then
       if(want_plot)then
          call modelcl%load(model_cl_file)
          ells = (/ (i, i=lmin_data, lmax_data) /)
          do i = 1, nmaps
             do j = i, nmaps
                ind = COOP_MATSYM_INDEX(nmaps, i, j)
                if(any(modelcl%cls(lmin_data:lmax_data, ind).ne.0.d0))then
                   call figure%open(trim(cl_output_root)//"_"//fieldnames(i:i)//fieldnames(j:j)//".txt")
                   call figure%init(xlabel = "$\ell$", ylabel = "$\ell(\ell+1)C_l^{"//fieldnames(i:i)//fieldnames(j:j)//"}/(2\pi)$")
                   call figure%plot(ells, modelcl%cls(lmin_data:lmax_data, ind)*ells*(ells+1.d0)/coop_2pi, color="red", linewidth=2.)
                   call figure%errorbars( x = templatepower%lb(1:num_ell_bins), y = mub(ind:(num_ell_bins-1)*numcls+ind:numcls)*templatepower%Cls(nint(templatepower%lb(1:num_ell_bins)), ind)* (templatepower%lb(1:num_ell_bins)*( templatepower%lb(1:num_ell_bins) + 1.d0)/coop_2pi), dy= abs(errorbars(ind:(num_ell_bins-1)*numcls+ind:numcls)) * templatepower%Cls(nint(templatepower%lb(1:num_ell_bins)), ind)* (templatepower%lb(1:num_ell_bins)*( templatepower%lb(1:num_ell_bins) + 1.d0)/coop_2pi),  color = "gray", center_color="black")
                   
                   call figure%close()
                endif
             enddo
          enddo
          call modelcl%free()
       endif
    endif

    do l = lmin_data, lmax_data
       do i = 1, numcls
          lambda = 0.d0
          do ib = templatepower%ibmin_l(l), templatepower%ibmax_l(l)
             lambda = lambda + mub(MUB_IND(ib, i))*templatepower%wl(ib, l)
          enddo
          templatepower%cls(l, i) = templatepower%cls(l, i) * lambda
       enddo
    enddo
    call templatepower%dump(trim(cl_output_root)//"_ITER"//COOP_STR_OF(iter)//".dat")

    do ind = 1, numcross
       call signalpower(ind)%free()
       call noisepower(ind)%free()
    enddo
    call templatepower%free()
    deallocate(index_mub_used, Fisher_used)
  end subroutine do_mub

  subroutine compare_pseudo_cls()
    COOP_STRING::noise_cl_root, data_cl_root, signal_cl_root
    COOP_INT::ich1, ich2, ind, i, j
    type(coop_asy)::figure
    type(coop_cls)::noise_cls, signal_cls, data_cls, model_cls
    COOP_REAL::ells(lmin_data:lmax_data)
    call coop_dictionary_lookup(settings, "signal_cl_root", signal_cl_root)
    call coop_dictionary_lookup(settings, "noise_cl_root", noise_cl_root)    
    call coop_dictionary_lookup(settings, "data_cl_root", data_cl_root)    
    ells = (/ (i, i=lmin_data, lmax_data) /)
    call model_cls%load(model_cl_file)
    call model_cls%filter(lmin = lmin_data, lmax = lmax_data)
    do ich1 = 1, num_channels
       do ich2 = ich1, num_channels
          call noise_cls%load(trim(noise_cl_root)//"_"//CHIND(ich1)//"_"//CHIND(ich2)//".fits")
          call signal_cls%load(trim(signal_cl_root)//"_"//CHIND(ich1)//"_"//CHIND(ich2)//".fits")
          call data_cls%load(trim(data_cl_root)//"_"//CHIND(ich1)//"_"//CHIND(ich2)//".fits")
          do i = 1, nmaps
             do j = i, nmaps
                ind = COOP_MATSYM_INDEX(nmaps, i, j)
                if(any(model_cls%cls(lmin_data:lmax_data, ind).ne.0.d0))then
                   call figure%open(trim(cl_output_root)//"_"//CHIND(ich1)//"_"//CHIND(ich2)//"_CMP_"//fieldnames(i:i)//fieldnames(j:j)//".txt")
                   call figure%init(xlabel = "$\ell$", ylabel = "$\ell(\ell+1)C_l^{"//fieldnames(i:i)//fieldnames(j:j)//"}/(2\pi)$")
                   call figure%plot(ells, (signal_cls%cls(lmin_data:lmax_data, ind)+noise_cls%cls(lmin_data:lmax_data, ind))*ells*(ells+1.d0)/coop_2pi, color="red", linewidth=1.5, legend = "Signal+Noise")
                   call figure%plot(ells, signal_cls%cls(lmin_data:lmax_data, ind)*ells*(ells+1.d0)/coop_2pi, color="black", linewidth=1., linetype="solid", legend = "Signal")
                   call figure%plot(ells, data_cls%cls(lmin_data:lmax_data, ind)*ells*(ells+1.d0)/coop_2pi, color="blue", linewidth=2., linetype="dotted", legend = "Data")
                   call figure%legend()
                   call figure%close()
                endif
             enddo
          enddo
       enddo
    enddo
  end subroutine compare_pseudo_cls


!!****************************************************************************

  subroutine load_qb(iter)
    COOP_INT::iter, ib
    COOP_REAL::junk
    call fp%open(trim(qb_output_root)//"_ITER"//COOP_STR_OF(iter)//".dat")
    do ib = 1, num_ell_bins
       read(fp%unit, *) junk, qb((ib-1)*numcross*numcls+1:ib*numcross*numcls)
    enddo
    call fp%close()
  end subroutine load_qb

  subroutine load_kernel(filename, lmin, lmax, kernel, weight)
    COOP_UNKNOWN_STRING::filename
    COOP_INT::lmin, lmax, lmin_saved, lmax_saved
    COOP_REAL::kernel(lmin:lmax, lmin:lmax, 4), weight
    COOP_REAL,dimension(:,:,:),allocatable::tmp
    type(coop_fits_file)::fp
    call fp%open_image(filename)
    call coop_dictionary_lookup(fp%header, "LMIN", lmin_saved)
    call coop_dictionary_lookup(fp%header, "LMAX", lmax_saved)
    call coop_dictionary_lookup(fp%header, "WEIGHT", weight)
    if(lmin .lt. lmin_saved .or. lmax .gt. lmax_saved)then
       stop "kernel ell range overflow; please recompute the kernel with sufficiently large lmax"
    endif
    allocate(tmp(lmin_saved:lmax_saved, lmin_saved:lmax_saved, 4))
    call fp%load_image_3d(tmp)
    kernel = tmp(lmin:lmax, lmin:lmax, :)
    call fp%close()
    deallocate(tmp)
  end subroutine load_kernel


  subroutine read_map(map, filename, maskfile)
    type(coop_healpix_maps)::map, mask
    COOP_UNKNOWN_STRING::filename
    COOP_UNKNOWN_STRING,optional::maskfile
    COOP_INT::i
    call map%read(filename, nmaps_wanted = nmaps, nested = .false.)
    map%map = max(min(map%map, map_maxval), map_minval)
    map%spin = spin
    if(present(maskfile))then
       call mask%read(maskfile, nmaps_wanted = 1, nested = .false.)
       do i = 1, nmaps
          map%map(:,i) = map%map(:,i)*mask%map(:,1)
       enddo
       call mask%free()
    endif
    call map%map2alm(lmax = lmax_data)
  end subroutine read_map


  subroutine get_Cbs(power, l, kernel, Cb_EE, Cb_BB)
    type(coop_binned_cls),intent(IN)::power
    COOP_INT,intent(IN)::l
    COOP_REAL,intent(IN)::kernel(power%lmin:power%lmax, power%lmin:power%lmax, 4)
    COOP_REAL,intent(OUT)::Cb_EE(numcls, num_ell_bins), Cb_BB(numcls, num_ell_bins)
    COOP_INT::ib, lp
    select case(trim(map_genre))
    case("I")
       Cb_EE = 0.d0
       do ib = 1, power%nb
          do lp = power%lmin_b(ib), power%lmax_b(ib)
             Cb_EE(1, ib) = Cb_EE(1, ib) + kernel(l, lp, coop_pseudoCl_kernel_index_TT)*power%Cls(lp, 1)*power%wl(ib, lp)
          enddo
       enddo
       Cb_BB = Cb_EE
    case("QU")
       Cb_EE = 0.d0
       Cb_BB = 0.d0
       do ib = 1, power%nb
          do lp = power%lmin_b(ib), power%lmax_b(ib)
             Cb_EE(coop_EB_index_EB, ib) = Cb_EE(coop_EB_index_EB, ib) + kernel(l, lp, coop_pseudoCl_kernel_index_EB) * power%Cls(lp, coop_EB_index_EB) * power%wl(ib, lp)
             Cb_EE(coop_EB_index_EE, ib) = Cb_EE(coop_EB_index_EE, ib) + (kernel(l, lp, coop_pseudoCl_kernel_index_EE_plus_BB) + kernel(l, lp, coop_pseudoCl_kernel_index_EE_minus_BB))/2.d0*power%Cls(lp, coop_EB_index_EE)*power%wl(ib, lp)
             Cb_EE(coop_EB_index_BB, ib) = Cb_EE(coop_EB_index_BB, ib) + (kernel(l, lp, coop_pseudoCl_kernel_index_EE_plus_BB) - kernel(l, lp, coop_pseudoCl_kernel_index_EE_minus_BB))/2.d0*power%Cls(lp, coop_EB_index_EE)*power%wl(ib, lp)
             Cb_BB(coop_EB_index_EE, ib) = Cb_BB(coop_EB_index_EE, ib) + (kernel(l, lp, coop_pseudoCl_kernel_index_EE_plus_BB) - kernel(l, lp, coop_pseudoCl_kernel_index_EE_minus_BB))/2.d0*power%Cls(lp, coop_EB_index_BB)*power%wl(ib, lp)
             Cb_BB(coop_EB_index_BB, ib) = Cb_BB(coop_EB_index_BB, ib) + (kernel(l, lp, coop_pseudoCl_kernel_index_EE_plus_BB) + kernel(l, lp, coop_pseudoCl_kernel_index_EE_minus_BB))/2.d0*power%Cls(lp, coop_EB_index_BB)*power%wl(ib, lp)
          enddo
       enddo
       Cb_BB(coop_EB_index_EB, :) =  Cb_EE(coop_EB_index_EB, :)
    case("IQU")
       Cb_EE = 0.d0
       Cb_BB = 0.d0
       do ib = 1, power%nb
          do lp = power%lmin_b(ib), power%lmax_b(ib)
             Cb_EE(coop_TEB_index_EB, ib) = Cb_EE(coop_TEB_index_EB, ib) + kernel(l, lp, coop_pseudoCl_kernel_index_EB) * power%Cls(lp, coop_TEB_index_EB) * power%wl(ib, lp)
             Cb_EE(coop_TEB_index_TT, ib) = Cb_EE(coop_TEB_index_TT, ib) + kernel(l, lp, coop_pseudoCl_kernel_index_TT) * power%Cls(lp, coop_TEB_index_TT) * power%wl(ib, lp)
             Cb_EE(coop_TEB_index_TE, ib) = Cb_EE(coop_TEB_index_TE, ib) + kernel(l, lp, coop_pseudoCl_kernel_index_TE) * power%Cls(lp, coop_TEB_index_TE) * power%wl(ib, lp)
             Cb_EE(coop_TEB_index_TB, ib) = Cb_EE(coop_TEB_index_TB, ib) + kernel(l, lp, coop_pseudoCl_kernel_index_TB) * power%Cls(lp, coop_TEB_index_TB) * power%wl(ib, lp)
             Cb_EE(coop_TEB_index_EE, ib) = Cb_EE(coop_TEB_index_EE, ib) + (kernel(l, lp, coop_pseudoCl_kernel_index_EE_plus_BB) + kernel(l, lp, coop_pseudoCl_kernel_index_EE_minus_BB))/2.d0*power%Cls(lp, coop_TEB_index_EE)*power%wl(ib, lp)
             Cb_EE(coop_TEB_index_BB, ib) = Cb_EE(coop_TEB_index_BB, ib) + (kernel(l, lp, coop_pseudoCl_kernel_index_EE_plus_BB) - kernel(l, lp, coop_pseudoCl_kernel_index_EE_minus_BB))/2.d0*power%Cls(lp, coop_TEB_index_EE)*power%wl(ib, lp)
             Cb_BB(coop_TEB_index_EE, ib) = Cb_BB(coop_TEB_index_EE, ib) + (kernel(l, lp, coop_pseudoCl_kernel_index_EE_plus_BB) - kernel(l, lp, coop_pseudoCl_kernel_index_EE_minus_BB))/2.d0*power%Cls(lp, coop_TEB_index_BB)*power%wl(ib, lp)
             Cb_BB(coop_TEB_index_BB, ib) = Cb_BB(coop_TEB_index_BB, ib) + (kernel(l, lp, coop_pseudoCl_kernel_index_EE_plus_BB) + kernel(l, lp, coop_pseudoCl_kernel_index_EE_minus_BB))/2.d0*power%Cls(lp, coop_TEB_index_BB)*power%wl(ib, lp)
          enddo
       enddo
       Cb_BB(coop_TEB_index_TT, :) =  Cb_EE(coop_TEB_index_TT, :)
       Cb_BB(coop_TEB_index_TE, :) =  Cb_EE(coop_TEB_index_TE, :)
       Cb_BB(coop_TEB_index_TB, :) =  Cb_EE(coop_TEB_index_TB, :)
       Cb_BB(coop_TEB_index_EB, :) =  Cb_EE(coop_TEB_index_EB, :)
    case default
       write(*,*) "map_genre = "//trim(map_genre)
       stop "XFASTER only supports map_genre = I/QU/IQU"
    end select
  end subroutine get_Cbs


  function sim_index_name(isim) result(str)
    COOP_STRING::str
    COOP_INT::isim
    if(sim_index_width .eq. 0)then
       str = "SIM"//COOP_STR_OF(isim)
    else
       str = coop_Ndigits(isim, sim_index_width)
    endif
  end function sim_index_name


end program Coop_Xfaster

