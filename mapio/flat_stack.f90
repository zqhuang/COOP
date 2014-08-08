program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  integer,parameter::lmin = 200
  integer,parameter::lmax = 2000
  integer,parameter::irepeat = 300
  COOP_REAL, parameter::reg_limit = 0.0005
  character(LEN=*),parameter::mapdir = "act/"
  character(LEN=*),parameter::mapid = "6"
  character(LEN=*),parameter::Ifile = mapdir//"I"//mapid//".fits"
  character(LEN=*),parameter::Qfile = mapdir//"Q"//mapid//".fits"
  character(LEN=*),parameter::Ufile = mapdir//"U"//mapid//".fits"
  character(LEN=*),parameter::Imaskfile = mapdir//"WI"//mapid//".fits"
  type(coop_fits_image_cea)::imap, umap, qmap, imask, cutmask
  type(coop_asy)::asy
  integer, parameter::n=300
  integer ix, iy, i, l
  type(coop_file) fp
  COOP_REAL, parameter::patchsize = 30.d0*coop_SI_arcmin
  COOP_UNKNOWN_STRING,parameter::output_dir = "ACTstacking/"

  COOP_REAL map(n, n), Cls(lmin:lmax)
  COOP_REAL, parameter::smooth_scale = coop_SI_arcmin * 1.5
  call imap%open(Ifile)
  call imask%open(Imaskfile)
  call qmap%open(Qfile)
  call umap%open(Ufile)

  where(imap%image .gt. 1000.)
     qmap%image = 0.
     umap%image = 0.
  end where

  where (imask%image .lt. 5000.)
     imap%image = imap%image * (imask%image/5000.)**8
     qmap%image = qmap%image * (imask%image/5000.)**8
     umap%image = umap%image * (imask%image/5000.)**8
  end where
  print*, maxval(imap%image), minval(imap%image)
  call imap%regularize(reg_limit)
  print*, maxval(imap%image), minval(imap%image)
  call imap%get_flatmap(smooth_scale)


  where(qmap%image**2+umap%image**2 .gt. 2.e5)
     qmap%image = 0.
     umap%image = 0.
  end where

  call qmap%get_flatmap(smooth_scale)
  call umap%get_flatmap(smooth_scale)
  call imap%get_QU(qmap, umap)
  call imap%smooth_flat(lmin = lmin, lmax = lmax)


  cutmask = imask
  where(imask%image .lt. 5000)
     cutmask%image = 0.
  elsewhere
     cutmask%image = 1.
  end where
  call cutmask%get_flatmap(smooth_scale)

!!$  qmap = imap
!!$  call qmap%get_QTUT()
  call imap%find_extrema(cutmask, "spots/I"//mapid//"_Tmax.txt", "Tmax", patchsize, irepeat)
!!$  call qmap%find_extrema(cutmask, "spots/I"//mapid//"_Tmax_QTUTOrient.txt", "Tmax_QTUTOrient", patchsize)
!!$  call qmap%find_extrema(cutmask, "spots/I"//mapid//"_PTmax.txt", "PTmax", patchsize)


!!$
!!$  call imap%stack2fig("spots/I"//mapid//"_Tmax.txt", "T", patchsize, output_dir//"stack_ACT_deep"//mapid//"_T_onTmax.txt", caption="$T$ on $T_{\max}$", label = "$T (\mu K)$", color_table = "Rainbow")
!!$  call imap%stack2fig("spots/I"//mapid//"_Tmax_QTUTOrient.txt", "T", patchsize, output_dir//"stack_ACT_deep"//mapid//"_T_onTmax_Oriented.txt", caption="$T$ on oriented $T_{\max}$", label = "$T (\mu K)$", color_table = "Rainbow")
!!$  call imap%stack2fig("spots/I"//mapid//"_Tmax_QTUTOrient.txt", "Q", patchsize, output_dir//"stack_ACT_deep"//mapid//"_Q_onTmax_Oriented.txt", caption="$Q$ on oriented $T_{\max}$", label = "$Q (\mu K)$", color_table = "Rainbow")
!!$  call imap%stack2fig("spots/I"//mapid//"_PTmax.txt", "Q", patchsize, output_dir//"stack_ACT_deep"//mapid//"_Q_onPTmax_Oriented.txt", caption="$Q$ on oriented $P_{T,\max}$", label = "$Q (\mu K)$", color_table = "Rainbow")
  call imap%stack2fig("spots/I"//mapid//"_Tmax.txt", "Qr", patchsize, output_dir//"stack_ACT_deep"//mapid//"_Qr_onTmax.txt", caption="$Q_r$ on $T_{\max}$", label = "$Q_r (\mu K)$", color_table = "Rainbow")


!!simulations

  call imap%simulate_flat(lmin = lmin,lmax = lmax, cls_file = "cls.dat")
  call imap%find_extrema(cutmask, "spots/simu_I"//mapid//"_Tmax.txt", "Tmax", patchsize, irepeat)

!!$  qmap = imap
!!$  call qmap%get_QTUT()
!!$  call qmap%find_extrema(cutmask, "spots/simu_I"//mapid//"_Tmax_QTUTOrient.txt", "Tmax_QTUTOrient", patchsize)
!!$
!!$  call qmap%find_extrema(cutmask, "spots/simu_I"//mapid//"_PTmax.txt", "PTmax", patchsize)


!!$  call imap%stack2fig("spots/simu_I"//mapid//"_Tmax.txt", "T", patchsize, output_dir//"stack_simu_deep"//mapid//"_T_onTmax.txt", caption="$T$ on $T_{\max}$", label = "$T (\mu K)$", color_table = "Rainbow")
!!$  call imap%stack2fig("spots/simu_I"//mapid//"_Tmax_QTUTOrient.txt", "T", patchsize, output_dir//"stack_simu_deep"//mapid//"_T_onTmax_Oriented.txt", caption="$T$ on oriented $T_{\max}$", label = "$T (\mu K)$", color_table = "Rainbow")
!!$
!!$  call imap%stack2fig("spots/simu_I"//mapid//"_Tmax_QTUTOrient.txt", "Q", patchsize, output_dir//"stack_simu_deep"//mapid//"_Q_onTmax_Oriented.txt", caption="$Q$ on oriented $T_{\max}$", label = "$Q (\mu K)$", color_table = "Rainbow")
!!$  call imap%stack2fig("spots/simu_I"//mapid//"_PTmax.txt", "Q", patchsize, output_dir//"stack_simu_deep"//mapid//"_Q_onPTmax_Oriented.txt" , caption="$Q$ on oriented $P_{T,\max}$", label = "$Q (\mu K)$", color_table = "Rainbow")
  call imap%stack2fig("spots/simu_I"//mapid//"_Tmax.txt", "Qr", patchsize, output_dir//"stack_simu_deep"//mapid//"_Qr_onTmax.txt", caption="$Q_r$ on $T_{\max}$", label = "$Q_r (\mu K)$", color_table = "Rainbow")


!!$  call imap%QU2EB()
!!$  call imap%stack2fig("spots/simu_I"//mapid//"_Tmax.txt", "E", patchsize, output_dir//"stack_simu_deep"//mapid//"_E_onTmax.txt", caption="$E$ on $T_{\max}$", label = "$E (\mu K)$", color_table = "Rainbow")


end program test
