program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  COOP_UNKNOWN_STRING,parameter::dir = "act16/"
  type(coop_fits_image_cea)::imap, qmap, umap, weight_II, weight_IQ, weight_IU, weight_QQ, weight_UU, weight_QU
  type(coop_fits_image_cea)::isum, qsum, usum, II_sum, IQ_sum, IU_sum, QQ_sum, UU_sum, QU_sum, mask
  COOP_INT:: i, istart, iend
  COOP_REAL::cov(3,3),  vec(3), det
  call coop_get_Input(1, istart)
  call coop_get_Input(2, iend)
  print*, istart, iend
  call mask%read(dir//"act_matcoadd_weight.fits")
  do i = istart, iend, iend-istart
     call imap%read(dir//"deep56_array_2_season2_iqu_c7v5_night_strict_nomoon_4way_set_"//COOP_STR_OF(i)//"_8Dec15_beams_srcsub_mapsub_wpoly_nobad_500_I.fits")
     call qmap%read(dir//"deep56_array_2_season2_iqu_c7v5_night_strict_nomoon_4way_set_"//COOP_STR_OF(i)//"_8Dec15_beams_srcsub_mapsub_wpoly_nobad_500_Q.fits")
     call umap%read(dir//"deep56_array_2_season2_iqu_c7v5_night_strict_nomoon_4way_set_"//COOP_STR_OF(i)//"_8Dec15_beams_srcsub_mapsub_wpoly_nobad_500_U.fits")
     call weight_II%read(dir//"deep56_array_2_season2_iqu_c7v5_night_strict_nomoon_4way_set_"//COOP_STR_OF(i)//"_8Dec15_beams_srcsub_mapsub_nobad_weights_I.fits")
     call weight_QQ%read(dir//"deep56_array_2_season2_iqu_c7v5_night_strict_nomoon_4way_set_"//COOP_STR_OF(i)//"_8Dec15_beams_srcsub_mapsub_nobad_weights_QQ.fits")
     call weight_UU%read(dir//"deep56_array_2_season2_iqu_c7v5_night_strict_nomoon_4way_set_"//COOP_STR_OF(i)//"_8Dec15_beams_srcsub_mapsub_nobad_weights_UU.fits")
     call weight_QU%read(dir//"deep56_array_2_season2_iqu_c7v5_night_strict_nomoon_4way_set_"//COOP_STR_OF(i)//"_8Dec15_beams_srcsub_mapsub_nobad_weights_QU.fits")
     call weight_IQ%read(dir//"deep56_array_2_season2_iqu_c7v5_night_strict_nomoon_4way_set_"//COOP_STR_OF(i)//"_8Dec15_beams_srcsub_mapsub_nobad_weights_Q.fits")
     call weight_IU%read(dir//"deep56_array_2_season2_iqu_c7v5_night_strict_nomoon_4way_set_"//COOP_STR_OF(i)//"_8Dec15_beams_srcsub_mapsub_nobad_weights_U.fits")
     call imap%regularize(0.003d0)
     call qmap%regularize(0.003d0)
     call umap%regularize(0.003d0)
     if(i.eq. istart)then
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
        II_sum%image = II_sum%image + weight_II%image
        QQ_sum%image  = QQ_sum%image + weight_QQ%image
        UU_sum%image =  UU_sum%image + weight_UU%image
        QU_sum%image = QU_sum%image + weight_QU%image
        IQ_sum%image = IQ_sum%image + weight_IQ%image
        IU_sum%image = IU_sum%image + weight_IU%image
        isum%image =  isum%image + imap%image * weight_II%image + qmap%image * weight_IQ%image + umap%image * weight_IU%image
        qsum%image = qsum%image + imap%image * weight_IQ%image + qmap%image * weight_QQ%image + umap%image * weight_QU%image
        usum%image =  usum%image + imap%image * weight_IU%image + qmap%image * weight_QU%image + umap%image * weight_UU%image
     endif
  enddo
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
     if(mask%image(i).gt.0.d0)then
        if(det .lt. 1.d3) stop "weird small det with large full det"
        cov = COOP_INV33(cov, det)
        vec = matmul(cov, (/ isum%image(i), qsum%image(i), usum%image(i) /) )
        imap%image(i) = vec(1)
        qmap%image(i) = vec(2)
        umap%image(i) = vec(3)
        weight_II%image(i) = det
     else
        weight_II%image(i) = 0.d0
        imap%image(i) = 0.d0
        qmap%image(i) = 0.d0
        umap%image(i) = 0.d0
     endif
  enddo
  write(*,"(A, F10.1,A)") "Fraction of sky used: ", 100.d0*count(weight_II%image .gt. 0.d0)/weight_II%npix, "%"
  call imap%write(dir//"act_matcoadd"//COOP_STR_OF(istart)//COOP_STR_OF(iend)//"_I.fits")
  call qmap%write(dir//"act_matcoadd"//COOP_STR_OF(istart)//COOP_STR_OF(iend)//"_Q.fits")
  call umap%write(dir//"act_matcoadd"//COOP_STR_OF(istart)//COOP_STR_OF(iend)//"_U.fits")
!  call weight_II%write(dir//"act_matcoadd_weight.fits")
end program test
