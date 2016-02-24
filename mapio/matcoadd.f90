program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  COOP_INT::nsets = 6
  COOP_UNKNOWN_STRING,parameter::dir = "act16/", out_prefix = "act_all"
  COOP_STRING::map_prefix, weight_prefix
  type(coop_fits_image_cea)::imap, qmap, umap, weight_II, weight_IQ, weight_IU, weight_QQ, weight_UU, weight_QU
  type(coop_fits_image_cea)::isum, qsum, usum, II_sum, IQ_sum, IU_sum, QQ_sum, UU_sum, QU_sum, mask
  COOP_INT:: i, pix, mpix, ix, iy, ixlow, ixup, iylow, iyup, ibase, mibase
  COOP_REAL::cov(3,3),  vec(3), det


  do i = 0, nsets - 1
     select case(trim(out_prefix))
     case("deep56_ar2")
        map_prefix = "deep56_array_2_season2_iqu_c7v5_night_strict_nomoon_4way_set_"//COOP_STR_OF(i)//"_8Dec15_beams_srcsub_mapsub_wpoly_nobad_500"
        weight_prefix = "deep56_array_2_season2_iqu_c7v5_night_strict_nomoon_4way_set_"//COOP_STR_OF(i)//"_8Dec15_beams_srcsub_mapsub_nobad_weights"
     case("deep56_ar1")
        map_prefix = "deep56_array_1_season2_iqu_c7v5_night_strict_nomoon_4way_set_"//COOP_STR_OF(i)//"_8Dec15_beams_srcsub_mapsub_wpoly_nobad_500"
        weight_prefix = "deep56_array_1_season2_iqu_c7v5_night_strict_nomoon_4way_set_"//COOP_STR_OF(i)//"_8Dec15_beams_srcsub_mapsub_nobad_weights"
     case("deep5_s1")
        map_prefix = "deep5_season1_4way_iqu_c7v5_"//COOP_STR_OF(i)//"_azfilt_strict_nomoon_srcpointing_noisecuts_srcsub2_2pass_wpoly_1000"
        weight_prefix = "deep5_season1_4way_iqu_c7v5_"//COOP_STR_OF(i)//"_azfilt_strict_nomoon_srcpointing_noisecuts_srcsub2_2pass_weights"
     case("deep6_s1")       
        map_prefix = "deep6_season1_4way_iqu_c7v5_"//COOP_STR_OF(i)//"_azfilt_strict_nomoon_srcpointing_noisecuts_srcsub2_2pass_wpoly_1000"
        weight_prefix = "deep6_season1_4way_iqu_c7v5_"//COOP_STR_OF(i)//"_azfilt_strict_nomoon_srcpointing_noisecuts_srcsub2_2pass_weights"
     case("deep6_ar1")
        map_prefix = "deep6_season2_4way_iqu_c7v5_arr1_"//COOP_STR_OF(i)//"_azfilt_strict_nomoon_srcpointing_noisecuts_srcsub2_2pass_wpoly_1000"
        weight_prefix = "deep6_season2_4way_iqu_c7v5_arr1_"//COOP_STR_OF(i)//"_azfilt_strict_nomoon_srcpointing_noisecuts_srcsub2_2pass_weights"
     case("deep6_ar2")
        map_prefix = "deep6_season2_4way_iqu_c7v5_arr2_"//COOP_STR_OF(i)//"_azfilt_strict_nomoon_srcpointing_noisecuts_srcsub2_2pass_wpoly_1000"
        weight_prefix = "deep6_season2_4way_iqu_c7v5_arr2_"//COOP_STR_OF(i)//"_azfilt_strict_nomoon_srcpointing_noisecuts_srcsub2_2pass_weights"
     case('deep56')
        map_prefix = "deep56_ar"//COOP_STR_OF(i+1)
        weight_prefix = "deep56_ar"//COOP_STR_OF(i+1)//"_weights"
     case("act_all")
        select case(i)
        case(0)
           map_prefix = "deep56_ar1"
           weight_prefix = "deep56_ar1_weights"
        case(1)
           map_prefix = "deep56_ar2"
           weight_prefix = "deep56_ar2_weights"
        case(2)
           map_prefix = "deep6_ar1"
           weight_prefix = "deep6_ar1_weights"
        case(3)
           map_prefix = "deep6_ar2"
           weight_prefix = "deep6_ar2_weights"
        case(4)
           map_prefix = "deep5_s1"
           weight_prefix = "deep5_s1_weights"
        case(5)
           map_prefix = "deep6_s1"
           weight_prefix = "deep6_s1_weights"
        end select
     case default
        stop "Unknown prefix"
     end select
     map_prefix = dir//trim(map_prefix)
     weight_prefix = dir//trim(weight_prefix)
     write(*,*) "coadding "//trim(map_prefix)
     call imap%read(trim(map_prefix)//"_I.fits")
     call qmap%read(trim(map_prefix)//"_Q.fits")
     call umap%read(trim(map_prefix)//"_U.fits")
     call weight_II%read(trim(weight_prefix)//"_I.fits")
     call weight_QQ%read(trim(weight_prefix)//"_QQ.fits")
     call weight_UU%read(trim(weight_prefix)//"_UU.fits")
     call weight_QU%read(trim(weight_prefix)//"_QU.fits")
     call weight_IQ%read(trim(weight_prefix)//"_Q.fits")
     call weight_IU%read(trim(weight_prefix)//"_U.fits")
     call imap%regularize(0.003d0)
     call qmap%regularize(0.003d0)
     call umap%regularize(0.003d0)
     if(i.eq.0)then
        mask = imap
        isum = imap
        qsum = qmap
        usum = umap
        II_sum = weight_II
        QQ_sum  = weight_QQ
        UU_sum = weight_UU
        QU_sum = weight_QU
        IQ_sum = weight_IQ
        IU_sum = weight_IU
        isum%image = imap%image * weight_II%image + qmap%image * weight_IQ%image + umap%image * weight_IU%image
        qsum%image = imap%image * weight_IQ%image + qmap%image * weight_QQ%image + umap%image * weight_QU%image
        usum%image = imap%image * weight_IU%image + qmap%image * weight_QU%image + umap%image * weight_UU%image
     else
        if(all(imap%nside .eq. isum%nside) .and. all(imap%icenter.eq.isum%icenter))then
           II_sum%image = II_sum%image + weight_II%image
           QQ_sum%image  = QQ_sum%image + weight_QQ%image
           UU_sum%image =  UU_sum%image + weight_UU%image
           QU_sum%image = QU_sum%image + weight_QU%image
           IQ_sum%image = IQ_sum%image + weight_IQ%image
           IU_sum%image = IU_sum%image + weight_IU%image
           isum%image =  isum%image + imap%image * weight_II%image + qmap%image * weight_IQ%image + umap%image * weight_IU%image
           qsum%image = qsum%image + imap%image * weight_IQ%image + qmap%image * weight_QQ%image + umap%image * weight_QU%image
           usum%image =  usum%image + imap%image * weight_IU%image + qmap%image * weight_QU%image + umap%image * weight_UU%image
        else
           write(*,*) "warning: you are coadding maps of different sizes"
           iylow = max(-imap%icenter(2), -isum%icenter(2))
           iyup = min(imap%nside(2)-imap%icenter(2), isum%nside(2)-isum%icenter(2))-1
           ixlow = max(-imap%icenter(1), -isum%icenter(1))
           ixup = min(imap%nside(1)-imap%icenter(1), isum%nside(1)-isum%icenter(1))-1
           !$omp parallel do private(pix, ix, iy, mpix, ibase, mibase)
           do iy = iylow, iyup
              ibase = (iy+imap%icenter(2))*imap%nside(1)+ imap%icenter(1)
              mibase = (iy+isum%icenter(2))*isum%nside(1) + isum%icenter(1)
              do ix =  ixlow, ixup
                 pix = ix  + ibase
                 mpix = ix + mibase
                 II_sum%image(mpix) = II_sum%image(mpix) + weight_II%image(pix)
                 QQ_sum%image(mpix) = QQ_sum%image(mpix) + weight_QQ%image(pix)
                 UU_sum%image(mpix) = UU_sum%image(mpix) + weight_UU%image(pix)
                 QU_sum%image(mpix) = QU_sum%image(mpix) + weight_QU%image(pix)
                 IQ_sum%image(mpix) = IQ_sum%image(mpix) + weight_IQ%image(pix)
                 IU_sum%image(mpix) = IU_sum%image(mpix) + weight_IU%image(pix)
                 isum%image(mpix) = isum%image(mpix) + imap%image(pix) * weight_II%image(pix) + qmap%image(pix) * weight_IQ%image(pix) + umap%image(pix) * weight_IU%image(pix)
                 qsum%image(mpix) = qsum%image(mpix) + imap%image(pix) * weight_IQ%image(pix) + qmap%image(pix) * weight_QQ%image(pix) + umap%image(pix) * weight_QU%image(pix)
                 usum%image(mpix) = usum%image(mpix) + imap%image(pix) * weight_IU%image(pix) + qmap%image(pix) * weight_QU%image(pix) + umap%image(pix) * weight_UU%image(pix)
              enddo
           enddo
           !$omp end parallel do
        endif
     endif
  enddo
  call imap%free()
  call qmap%free()
  call umap%free()
  call weight_II%free()
  call weight_IQ%free()
  call weight_IU%free()
  call weight_QQ%free()
  call weight_UU%free()
  call weight_QU%free()
  if(nsets .gt. 1)then
     do i = 0, isum%npix-1
        cov(1,1) = II_sum%image(i)
        cov(1,2) = IQ_sum%image(i)
        cov(1,3) = IU_sum%image(i)
        cov(2,2) = QQ_sum%image(i)
        cov(2,3) = QU_sum%image(i)
        cov(3,3) = UU_sum%image(i)
        cov(2,1) = cov(1,2)
        cov(3,1) = cov(1,3)
        cov(3,2) = cov(2,3)
        det = COOP_DET33(cov)
        if(det .gt. 1.d11)then
           cov = COOP_INV33(cov, det)
           vec = matmul(cov, (/ isum%image(i), qsum%image(i), usum%image(i) /) )
           isum%image(i) = vec(1)
           qsum%image(i) = vec(2)
           usum%image(i) = vec(3)
           mask%image(i) = det
        else
           mask%image(i) = 0.d0
           isum%image(i) = 0.d0
           qsum%image(i) = 0.d0
           usum%image(i) = 0.d0
           II_sum%image(i) = 0.d0
           IQ_sum%image(i) = 0.d0
           IU_sum%image(i) = 0.d0
           QQ_sum%image(i) = 0.d0
           UU_sum%image(i) = 0.d0
           QU_sum%image(i) = 0.d0
        endif
     enddo
  else
     do i = 0, isum%npix-1
        cov(1,1) = II_sum%image(i)
        cov(1,2) = IQ_sum%image(i)
        cov(1,3) = IU_sum%image(i)
        cov(2,2) = QQ_sum%image(i)
        cov(2,3) = QU_sum%image(i)
        cov(3,3) = UU_sum%image(i)
        cov(2,1) = cov(1,2)
        cov(3,1) = cov(1,3)
        cov(3,2) = cov(2,3)
        det = COOP_DET33(cov)
        if(det .gt. 1.d11)then
           mask%image(i) = det
        else
           mask%image(i) = 0.d0
        endif
     enddo
  endif
  call isum%write(dir//trim(out_prefix)//"_I.fits")
  call qsum%write(dir//trim(out_prefix)//"_Q.fits")
  call usum%write(dir//trim(out_prefix)//"_U.fits")

  call II_sum%write(dir//trim(out_prefix)//"_weights_I.fits")
  call QQ_sum%write(dir//trim(out_prefix)//"_weights_QQ.fits")
  call UU_sum%write(dir//trim(out_prefix)//"_weights_UU.fits")
  call QU_sum%write(dir//trim(out_prefix)//"_weights_QU.fits")
  call IQ_sum%write(dir//trim(out_prefix)//"_weights_Q.fits")
  call IU_sum%write(dir//trim(out_prefix)//"_weights_U.fits")

  call mask%write(dir//trim(out_prefix)//"_mask.fits")
end program test
