program coop_Xfaster
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use coop_fitsio_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  implicit none
#include "constants.h"
#define QB_IND(ib, icross, icl) (((ib-1)*numcross + icross - 1 )*numcls + icl)

!!this is coop version of xfaster, written from scratch. 
  COOP_STRING::inifile, qb_input_file, qb_output_root
  COOP_STRING::action, map_genre, map_unit
  COOP_INT::ib, num_iterations, iter
  type(coop_dictionary)::settings
  type(coop_file)::fp
  COOP_INT::num_ell_bins, num_channels, lmin_data, lmax_data, nmaps, numcls, numcross, num_qbs, matdim
  COOP_INT,dimension(:),allocatable::spin
  COOP_REAL,dimension(:),allocatable::qb

  if(iargc().lt.1)then
     write(*,*) "Syntax:"
     write(*,*) "./XFASTER inifile"
     stop
  endif
  inifile = trim(coop_InputArgs(1))
  call coop_load_dictionary(inifile, settings)
  call coop_dictionary_lookup(settings, 'action', action)
  call coop_dictionary_lookup(settings, 'map_genre', map_genre, 'IQU')
  call coop_dictionary_lookup(settings, 'qb_input_file', qb_input_file, "")
  call coop_dictionary_lookup(settings, 'qb_output_root', qb_output_root)
  call coop_dictionary_lookup(settings, 'num_iterations', num_iterations)

  nmaps = len_trim(map_genre)
  numcls = nmaps*(nmaps+1)/2
  
  allocate(spin(nmaps))
  select case(trim(map_genre))
  case("I")
    spin = 0
 case("QU")
    spin = (/ 2, 2 /)
 case("IQU")
    spin = (/ 0, 2, 2 /)
 case default
    write(*,*) "map_genre = "//trim(map_genre)
    stop "XFASTER only supports map_genre = I/QU/IQU"
 end select
     
  call coop_dictionary_lookup(settings, 'map_unit', map_unit, 'muK')

  call coop_dictionary_lookup(settings, 'num_channels', num_channels)
  numcross = num_channels*(num_channels+1)/2
  matdim = num_channels*nmaps

  call coop_dictionary_lookup(settings, "num_ell_bins", num_ell_bins)
  num_qbs = numcross*numcls*num_ell_bins

  allocate(qb(num_qbs))
  if(trim(qb_input_file).eq."")then
     qb = 1.d0
  else
     write(*,*) "loading "//trim(qb_input_file)
     call fp%open(qb_input_file)
     do ib = 1, num_ell_bins
        read(fp%unit, "("//COOP_STR_OF(numcross*numcls)//"E16.7)") qb((ib-1)*numcross*numcls+1:ib*numcross*numcls)
     enddo
     call fp%close()
  endif

  call coop_dictionary_lookup(settings, 'lmin_data', lmin_data, 2)
  call coop_dictionary_lookup(settings, 'lmax_data', lmax_data)
  coop_healpix_want_cls = .true.
  if(num_channels .lt. 1)then
     write(*,*) "number of channels = ", num_channels
     stop "XFASTER only supports num_channels >=1"
  endif
  select case(trim(action))
  case("DO_KERNEL")
     call do_kernel()
  case("DO_NOISE")
     call do_noise()
  case("DO_TRANSFER")
     call do_transfer()
  case("DO_PSEUDO_CL")
     call do_pseudo_cl()
  case("DO_CL")
     do iter = 1, num_iterations
        call do_cl(iter)
     enddo
  case("DO_ALL")
     call do_kernel()
     call do_noise()
     call do_transfer()
     call do_pseudo_cl()
     do iter = 1, num_iterations
        call do_cl(iter)
     enddo
  case default
     write(*,*) "action  = "//trim(action)
     stop "Unknown action"
  end select

contains

  subroutine do_kernel()
    COOP_STRING::kernel_root, mask_root
    COOP_INT::ich1, ich2, lmax_mask, lmax_kernel
    type(coop_healpix_maps)::mask1, mask2
    type(coop_dictionary)::header
    COOP_REAL,dimension(:,:,:),allocatable::kernel
    write(*,*) "***************************************************"
    write(*,*) "Doing mask kernels."
    call coop_dictionary_lookup(settings, "kernel_root", kernel_root)
    call coop_dictionary_lookup(settings, "mask_root", mask_root)
    call coop_dictionary_lookup(settings, "lmax_mask", lmax_mask, 300)
    call coop_dictionary_lookup(settings, "lmax_kernel", lmax_kernel, lmax_data)
    allocate(kernel(0:lmax_kernel, 0:lmax_kernel, 4))
    call header%insert("LMIN", "0")
    call header%insert("LMAX", COOP_STR_OF(lmax_kernel))
    do ich1=1, num_channels
       call mask1%read(trim(mask_root)//"_"//COOP_STR_OF(ich1)//".fits", nmaps_wanted=1, nested = .false.)
       call mask1%map2alm(lmax = lmax_mask)
       call header%insert("WEIGHT", COOP_STR_OF( sum(mask1%map(:,1)**2)**2/dble(mask1%npix)/sum(mask1%map(:,1)**4)), overwrite = .true.)
       call coop_pseudoCl_get_kernel(lmax_mask = lmax_mask, Cl_mask = dble(mask1%cl(0:lmax_mask, 1)), lmin = 0 , lmax= lmax_kernel, kernel = kernel)
       call coop_fits_file_write_image_3d(filename = trim(kernel_root)//"_"//COOP_STR_OF(ich1)//"_"//COOP_STR_OF(ich1)//".fits", image = kernel, header=header)
       write(*,*) "kernel file "//trim(kernel_root)//"_"//COOP_STR_OF(ich1)//"_"//COOP_STR_OF(ich1)//".fits is produced."
       do ich2 = ich1+1, num_channels
          call mask2%read(trim(mask_root)//"_"//COOP_STR_OF(ich2)//".fits", nmaps_wanted=1, nested = .false.)
          call mask2%map2alm(lmax = lmax_mask)
          call mask2%get_cls(mask1)
          call coop_pseudoCl_get_kernel(lmax_mask = lmax_mask, Cl_mask = dble(mask2%cl(0:lmax_mask, 1)), lmin = 0 , lmax= lmax_kernel, kernel = kernel)
          call header%insert("WEIGHT", COOP_STR_OF(sum(mask1%map(:,1)*mask2%map(:,1))**2/mask1%npix/sum((mask1%map(:,1)*mask2%map(:,1))**2)), overwrite = .true.)
          call coop_fits_file_write_image_3d(filename = trim(kernel_root)//"_"//COOP_STR_OF(ich1)//"_"//COOP_STR_OF(ich2)//".fits", image = kernel, header = header)
          write(*,*) "kernel file "//trim(kernel_root)//"_"//COOP_STR_OF(ich1)//"_"//COOP_STR_OF(ich2)//".fits is produced."
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
    write(*,*) "***************************************************"
    write(*,*) "doing noise power spectrum N_l from "//COOP_STR_OF(num_noise_sims)//" simulations"
    do ich1 = 1, num_channels
       do ich2 = ich1, num_channels
          write(*,"(A$)") "channel "//COOP_STR_OF(ich1)//" x channel "//COOP_STR_OF(ich2)
          noisepower%cls = 0.d0
          do isim = 1, num_noise_sims
             write(*,"(A$)") "."
             call read_map(noise1, trim(noise_sim_root)//"_SIM"//COOP_STR_OF(isim)//"_"//COOP_STR_OF(ich1)//".fits", trim(mask_root)//"_"//COOP_STR_OF(ich1)//".fits")
             if(ich2.ne.ich1)then
                call read_map(noise2, trim(noise_sim_root)//"_SIM"//COOP_STR_OF(isim)//"_"//COOP_STR_OF(ich2)//".fits", trim(mask_root)//"_"//COOP_STR_OF(ich2)//".fits")
                call noise1%get_cls(noise2)
             endif
             noisepower%cls = noisepower%cls + noise1%cl(lmin_data:lmax_data, :)
          enddo
          noisepower%cls = noisepower%cls/num_noise_sims
          call noisepower%smooth(5)
          call noisepower%dump(trim(noise_cl_root)//"_"//COOP_STR_OF(ich1)//"_"//COOP_STR_OF(ich2)//".fits")
          write(*,*)
       enddo
    enddo
    call noise1%free()
    call noise2%free()
    call noisepower%free()
  end subroutine do_noise

  subroutine do_transfer()
    COOP_STRING::model_cl_root, signal_sim_root, transfer_root
    type(coop_cls)::modelpower, signalpower
    COOP_INT::ich1, ich2, num_signal_sims, isim, ib
    type(coop_healpix_maps)::sig1, sig2
    call coop_dictionary_lookup(settings, "signal_sim_root", signal_sim_root)
    call coop_dictionary_lookup(settings, "model_cl_root", model_cl_root)
    call coop_dictionary_lookup(settings, "transfer_root", transfer_root)
    call coop_dictionary_lookup(settings, "num_signal_sims", num_signal_sims)
    write(*,*) "***************************************************"
    write(*,*) "doing transfer F_l from "//COOP_STR_OF(num_signal_sims)//" simulations"
    call signalpower%init(lmin = lmin_data, lmax = lmax_data, genre = map_genre, unit = map_unit, spin = spin)
    do ich1 = 1, num_channels
       do ich2 = ich1, num_channels
          write(*,"(A$)") "channel "//COOP_STR_OF(ich1)//" x channel "//COOP_STR_OF(ich2)
          call modelpower%load(trim(model_cl_root)//"_"//COOP_STR_OF(ich1)//"_"//COOP_STR_OF(ich2)//".fits")
          call modelpower%filter(lmin = lmin_data, lmax = lmax_data)
          signalpower%cls = 0.d0
          do isim = 1, num_signal_sims
             write(*,"(A$)") "."
             call read_map(sig1, trim(signal_sim_root)//"_SIM"//COOP_STR_OF(isim)//"_"//COOP_STR_OF(ich1)//".fits")
             if(ich2.ne.ich1)then
                call read_map(sig2, trim(signal_sim_root)//"_SIM"//COOP_STR_OF(isim)//"_"//COOP_STR_OF(ich2)//".fits")
                call sig1%get_cls(sig2)
             endif
             signalpower%cls = signalpower%cls + sig1%cl(lmin_data:lmax_data, :)
          enddo
          write(*,*)
          where(modelpower%cls .ne. 0)
             signalpower%cls = signalpower%cls/num_signal_sims/modelpower%cls
          elsewhere
             signalpower%cls = 0.d0
          end where
          call signalpower%smooth(5)
          call signalpower%dump(trim(transfer_root)//"_"//COOP_STR_OF(ich1)//"_"//COOP_STR_OF(ich2)//".fits")

       enddo
    enddo
    call sig1%free()
    call sig2%free()
    call signalpower%free()
    call modelpower%free()
  end subroutine do_transfer

  subroutine do_pseudo_cl()
    COOP_STRING::data_map_root, data_cl_root, mask_root
    type(coop_healpix_maps)::map1, map2
    COOP_INT::ich1, ich2
    type(coop_cls)::datapower
    write(*,*) "***************************************************"
    write(*,"(A$)") "Doing pseudo Cls"
    call coop_dictionary_lookup(settings, "data_map_root", data_map_root)
    call coop_dictionary_lookup(settings, "data_cl_root", data_cl_root)
    call coop_dictionary_lookup(settings, "mask_root", mask_root)
    do ich1 =1, num_channels
       call read_map(map1, trim(data_map_root)//"_"//COOP_STR_OF(ich1)//".fits", trim(mask_root)//"_"//COOP_STR_OF(ich1)//".fits")
       call datapower%init(lmin = lmin_data, lmax = lmax_data, genre = map_genre, unit = map_unit, spin = spin, cls = dble(map1%cl(lmin_data:lmax_data, :)))
       call datapower%dump(trim(data_cl_root)//"_"//COOP_STR_OF(ich1)//"_"//COOP_STR_OF(ich1)//".fits")
       write(*,"(A$)") "."
       do ich2 = ich1+1, num_channels
          call read_map(map2, trim(data_map_root)//"_"//COOP_STR_OF(ich2)//".fits", trim(mask_root)//"_"//COOP_STR_OF(ich2)//".fits")
          call map2%get_cls(map1)
          call datapower%init(lmin = lmin_data, lmax = lmax_data, genre = map_genre, unit = map_unit, spin = spin, cls = dble(map2%cl(lmin_data:lmax_data, :)))
          call datapower%dump(trim(data_cl_root)//"_"//COOP_STR_OF(ich1)//"_"//COOP_STR_OF(ich2)//".fits")
          write(*,"(A$)") "."
       enddo
    enddo
    write(*,*)
    call map1%free()
    call map2%free()
    call datapower%free()
  end subroutine do_pseudo_cl

  subroutine do_cl(iter)
    COOP_STRING::data_cl_root, noise_cl_root, template_cl_root, kernel_root, transfer_root
    COOP_INT:: l, iter
    type(coop_binned_cls),dimension(numcross)::datapower, noisepower, templatepower, trans
    COOP_REAL::kernel(lmin_data:lmax_data, lmin_data:lmax_data, 4,numcross)
    COOP_REAL::weight(numcross)
    COOP_INT::i, j, ind, ii, jj, ib, ich1, ich2, icl, jb, jch1, jch2, jcl, iqb
    COOP_REAL::mat(matdim, matdim), dSdq(matdim, matdim, num_qbs), wmat(matdim, matdim), datamat(matdim, matdim), invDdSdq(matdim, matdim, num_qbs)
    logical::nonzero(num_qbs)
    COOP_REAL::Fisher(num_qbs, num_qbs)
    COOP_REAL::vecb(num_qbs)
    COOP_REAL::Cb_EE(numcls, num_ell_bins), Cb_BB(numcls, num_ell_bins)
    type(coop_file)::fp
    write(*,*) "***************************************************"
    write(*,*) "Iterating the maximum-likelihood Cls"
    call coop_dictionary_lookup(settings, "kernel_root", kernel_root)
    call coop_dictionary_lookup(settings, "transfer_root", transfer_root)
    call coop_dictionary_lookup(settings, "data_cl_root", data_cl_root)
    call coop_dictionary_lookup(settings, "noise_cl_root", noise_cl_root)
    call coop_dictionary_lookup(settings, "template_cl_root", template_cl_root)
    do i = 1, num_channels
       do j = i, num_channels
          ind = COOP_MATSYM_INDEX(num_channels, i,j)
          call load_kernel(trim(kernel_root)//"_"//COOP_STR_OF(i)//"_"//COOP_STR_OF(j)//".fits", lmin_data, lmax_data, kernel(:,:,:,ind), weight(ind))
          call datapower(ind)%load(trim(data_cl_root)//"_"//COOP_STR_OF(i)//"_"//COOP_STR_OF(j)//".fits")
          call datapower(ind)%filter(lmin = lmin_data, lmax = lmax_data)
          call noisepower(ind)%load(trim(noise_cl_root)//"_"//COOP_STR_OF(i)//"_"//COOP_STR_OF(j)//".fits")
          call noisepower(ind)%filter(lmin = lmin_data, lmax = lmax_data)
          call trans(ind)%load(trim(transfer_root)//"_"//COOP_STR_OF(i)//"_"//COOP_STR_OF(j)//".fits")
          call trans(ind)%filter(lmin = lmin_data, lmax = lmax_data)
          call templatepower(ind)%load(trim(template_cl_root)//"_"//COOP_STR_OF(i)//"_"//COOP_STR_OF(j)//".fits")
          call templatepower(ind)%filter(lmin = lmin_data, lmax = lmax_data)
          templatepower(ind)%cls = templatepower(ind)%cls * trans(ind)%cls
          call templatepower(ind)%alloc(num_ell_bins)
       enddo
    enddo
    vecb = 0.d0
    Fisher = 0.d0
    do l = lmin_data, lmax_data
       dSdq = 0.d0
       do ich1 = 1, num_channels
          do ich2 = ich1, num_channels
             ind = COOP_MATSYM_INDEX(num_channels, ich1, ich2)
             call get_Cbs(templatepower(ind), l, kernel(:,:,:,ind), Cb_EE, Cb_BB)
             select case(trim(map_genre))
             case("I")
                do ib = 1, num_ell_bins
                   dSdq(ich1, ich2, QB_IND(ib, ind, 1)) = Cb_EE(1, ib)
                enddo
             case("QU")
                do ib = 1, num_ell_bins
                   dSdq((ich1-1)*nmaps+coop_EB_index_E, (ich2-1)*nmaps+coop_EB_index_E, QB_IND(ib, ind, coop_EB_index_EE)) = Cb_EE(coop_EB_index_EE, ib)
                   dSdq((ich1-1)*nmaps+coop_EB_index_B, (ich2-1)*nmaps+coop_EB_index_B, QB_IND(ib, ind, coop_EB_index_EE)) = Cb_EE(coop_EB_index_BB, ib)

                   dSdq((ich1-1)*nmaps+coop_EB_index_E, (ich2-1)*nmaps+coop_EB_index_E, QB_IND(ib, ind, coop_EB_index_BB)) = Cb_BB(coop_EB_index_EE, ib)
                   dSdq((ich1-1)*nmaps+coop_EB_index_B, (ich2-1)*nmaps+coop_EB_index_B, QB_IND(ib, ind, coop_EB_index_BB)) = Cb_BB(coop_EB_index_BB, ib)

                   dSdq((ich1-1)*nmaps+coop_EB_index_E, (ich2-1)*nmaps+coop_EB_index_B, QB_IND(ib, ind, coop_EB_index_EB)) = Cb_EE(coop_EB_index_EB, ib)
                   dSdq((ich1-1)*nmaps+coop_EB_index_B, (ich2-1)*nmaps+coop_EB_index_E, QB_IND(ib, ind, coop_EB_index_EB)) = Cb_EE(coop_EB_index_EB, ib)
                   
                enddo
             case("IQU")
                do ib = 1, num_ell_bins
                   dSdq((ich1-1)*nmaps+coop_TEB_index_E, (ich2-1)*nmaps+coop_TEB_index_E, QB_IND(ib, ind, coop_TEB_index_EE)) = Cb_EE(coop_TEB_index_EE, ib)
                   dSdq((ich1-1)*nmaps+coop_TEB_index_B, (ich2-1)*nmaps+coop_TEB_index_B, QB_IND(ib, ind, coop_TEB_index_EE)) = Cb_EE(coop_TEB_index_BB, ib)

                   dSdq((ich1-1)*nmaps+coop_TEB_index_E, (ich2-1)*nmaps+coop_TEB_index_E, QB_IND(ib, ind, coop_TEB_index_BB)) = Cb_BB(coop_TEB_index_EE, ib)
                   dSdq((ich1-1)*nmaps+coop_TEB_index_B, (ich2-1)*nmaps+coop_TEB_index_B, QB_IND(ib, ind, coop_TEB_index_BB)) = Cb_BB(coop_TEB_index_BB, ib)

                   dSdq((ich1-1)*nmaps+coop_TEB_index_E, (ich2-1)*nmaps+coop_TEB_index_B, QB_IND(ib, ind, coop_TEB_index_EB)) = Cb_EE(coop_TEB_index_EB, ib)
                   dSdq((ich1-1)*nmaps+coop_TEB_index_B, (ich2-1)*nmaps+coop_TEB_index_E, QB_IND(ib, ind, coop_TEB_index_EB)) = Cb_EE(coop_TEB_index_EB, ib)

                   dSdq((ich1-1)*nmaps+coop_TEB_index_T, (ich2-1)*nmaps+coop_TEB_index_T, QB_IND(ib, ind, coop_TEB_index_TT)) = Cb_EE(coop_TEB_index_TT, ib)
                   dSdq((ich1-1)*nmaps+coop_TEB_index_T, (ich2-1)*nmaps+coop_TEB_index_E, QB_IND(ib, ind, coop_TEB_index_TE)) = Cb_EE(coop_TEB_index_TE, ib)
                   dSdq((ich1-1)*nmaps+coop_TEB_index_T, (ich2-1)*nmaps+coop_TEB_index_B, QB_IND(ib, ind, coop_TEB_index_TB)) = Cb_EE(coop_TEB_index_TB, ib)

                   dSdq((ich1-1)*nmaps+coop_TEB_index_E, (ich2-1)*nmaps+coop_TEB_index_T, QB_IND(ib, ind, coop_TEB_index_TE)) = Cb_EE(coop_TEB_index_TE, ib)
                   dSdq((ich1-1)*nmaps+coop_TEB_index_B, (ich2-1)*nmaps+coop_TEB_index_T, QB_IND(ib, ind, coop_TEB_index_TB)) = Cb_EE(coop_TEB_index_TB, ib)
                enddo
             end select

             do i = 1, nmaps
                do j = i, nmaps
                   mat((ich1-1)*nmaps + i, (ich2-1)*nmaps+j) = noisepower(ind)%cls(l, COOP_MATSYM_INDEX(nmaps, i, j))
                   datamat((ich1-1)*nmaps + i, (ich2-1)*nmaps+j) = datapower(ind)%cls(l, COOP_MATSYM_INDEX(nmaps, i, j)) - noisepower(ind)%cls(l, COOP_MATSYM_INDEX(nmaps, i, j))
                   if(i.ne.j)then
                      mat((ich1-1)*nmaps + j, (ich2-1)*nmaps+i) = mat((ich1-1)*nmaps + i, (ich2-1)*nmaps + j)
                      datamat((ich1-1)*nmaps + j, (ich2-1)*nmaps+i) = datamat((ich1-1)*nmaps + i, (ich2-1)*nmaps + j)
                   endif

                enddo
             enddo
             if(ich1 .ne. ich2)then
                dSdq((ich2-1)*nmaps+1:ich2*nmaps, (ich1-1)*nmaps+1:ich1*nmaps, :) = dSdq((ich1-1)*nmaps+1:ich1*nmaps, (ich2-1)*nmaps+1:ich2*nmaps, :)
                mat((ich2-1)*nmaps+1:ich2*nmaps, (ich1-1)*nmaps+1:ich1*nmaps) = mat((ich1-1)*nmaps+1:ich1*nmaps, (ich2-1)*nmaps+1:ich2*nmaps)
                datamat((ich2-1)*nmaps+1:ich2*nmaps, (ich1-1)*nmaps+1:ich1*nmaps) = datamat((ich1-1)*nmaps+1:ich1*nmaps, (ich2-1)*nmaps+1:ich2*nmaps)
             endif
          enddo
       enddo
       do ib = 1, num_qbs
          if(all(dSdq(:,:,ib).eq.0.d0))then
             nonzero(ib) = .false.
          else
             mat = mat + dSdq(:,:, ib)*qb(ib)
             nonzero(ib) = .true.
          endif
       enddo
       call coop_sympos_inverse(matdim, matdim, mat)
       wmat = mat
!!$       do ich1 = 1, num_channels
!!$          do ich2 = 1, num_channels
!!$             wmat((ich1-1)*nmaps+1:ich1*nmaps, (ich2-1)*nmaps+1:ich2*nmaps) = mat((ich1-1)*nmaps+1:ich1*nmaps, (ich2-1)*nmaps+1:ich2*nmaps) *weight(COOP_MATSYM_INDEX(num_channels, ich1, ich2))
!!$          enddo
!!$       enddo
       do ib = 1, num_qbs
          if(nonzero(ib))then
             invDdSdq(:, :, ib) = matmul(wmat, matmul(dSdq(:,:, ib), mat))
          endif
       enddo
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
    enddo
    write(*,*)
    do ib = 1, num_qbs
       do jb = ib+1, num_qbs
          Fisher(jb, ib)  = Fisher(ib, jb)
       enddo
    enddo
    call coop_sympos_inverse(num_qbs, num_qbs, Fisher)
    qb = matmul(Fisher, vecb)
    call fp%open(trim(qb_output_root)//"_ITER"//COOP_STR_OF(iter)//".dat")
    do ib = 1, num_ell_bins
       write(fp%unit, "(F10.2, "//COOP_STR_OF(numcross*numcls)//"E16.7)") templatepower(1)%lb(ib), qb((ib-1)*numcross*numcls+1:ib*numcross*numcls)
    enddo
    call fp%close()

    do ind = 1, numcross
       call datapower(ind)%free()
       call noisepower(ind)%free()
       call templatepower(ind)%free()
       call trans(ind)%free()
    enddo
  end subroutine do_cl

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
    end select
  end subroutine get_Cbs



end program Coop_Xfaster

