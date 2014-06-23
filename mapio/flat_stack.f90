program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  integer,parameter::lmin = 200
  integer,parameter::lmax = 2500
  character(LEN=*),parameter::mapdir = "act/"
  character(LEN=*),parameter::Ifile = mapdir//"I6.fits"
  character(LEN=*),parameter::Qfile = mapdir//"Q6.fits"
  character(LEN=*),parameter::Ufile = mapdir//"U6.fits"
  character(LEN=*),parameter::Imaskfile = mapdir//"WI6.fits"
  type(coop_fits_image_cea)::imap, umap, qmap, imask
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
  call imap%regularize(0.005d0)
  where (imask%image .lt. 5000.)
     imap%image = imap%image * (imask%image/5000.)**8
  end where
  call imap%get_flatmap(smooth_scale)


  call qmap%open(Qfile)
  call qmap%regularize(0.002d0)
  where (imask%image .lt. 5000.)
     qmap%image = qmap%image * (imask%image/5000.)**8
  end where
  call qmap%get_flatmap(smooth_scale)

  call umap%open(Ufile)
  call umap%regularize(0.002d0)
  where (imask%image .lt. 5000.)
     umap%image = umap%image * (imask%image/5000.)**8
  end where
  call umap%get_flatmap(smooth_scale)

  call imap%get_QU(qmap, umap)
  call imap%smooth_flat(lmin = lmin, lmax = lmax)



  where(imask%image .lt. 5000)
     imask%image = 0.
  elsewhere
     imask%image = 1.
  end where
  call imask%get_flatmap(smooth_scale)

  qmap = imap
  call qmap%get_QTUT()
  call qmap%find_extrema(imask, "spots/I6_Tmax.txt", "Tmax", patchsize)
  call qmap%find_extrema(imask, "spots/I6_Tmax_QTUTOrient.txt", "Tmax_QTUTOrient", patchsize)
  call qmap%find_extrema(imask, "spots/I6_PTmax.txt", "PTmax", patchsize)



  call imap%stack2fig("spots/I6_Tmax.txt", "T", patchsize, output_dir//"stack_ACT_deep6_T_onTmax.txt", caption="$T$ on $T_{\max}$", label = "$T (\mu K)$", color_table = "Rainbow")
  call imap%stack2fig("spots/I6_Tmax_QTUTOrient.txt", "T", patchsize, output_dir//"stack_ACT_deep6_T_onTmax_Oriented.txt", caption="$T$ on oriented $T_{\max}$", label = "$T (\mu K)$", color_table = "Rainbow")
  call imap%stack2fig("spots/I6_Tmax_QTUTOrient.txt", "Q", patchsize, output_dir//"stack_ACT_deep6_Q_onTmax_Oriented.txt", caption="$Q$ on oriented $T_{\max}$", label = "$Q (\mu K)$", color_table = "Rainbow")
  call imap%stack2fig("spots/I6_PTmax.txt", "Q", patchsize, output_dir//"stack_ACT_deep6_Q_onPTmax_Oriented.txt", caption="$Q$ on oriented $P_{T,\max}$", label = "$Q (\mu K)$", color_table = "Rainbow")
  call imap%stack2fig("spots/I6_Tmax.txt", "Qr", patchsize, output_dir//"stack_ACT_deep6_Qr_onTmax.txt", caption="$Q_r$ on $T_{\max}$", label = "$Q_r (\mu K)$", color_table = "Rainbow")


!!simulations

  call imap%simulate_flat(lmin = lmin,lmax = lmax, cls_file = "cls.txt")
  call imap%find_extrema(imask, "spots/simu_I6_Tmax.txt", "Tmax", patchsize)

  qmap = imap
  call qmap%get_QTUT()
  call qmap%find_extrema(imask, "spots/simu_I6_Tmax_QTUTOrient.txt", "Tmax_QTUTOrient", patchsize)

  call qmap%find_extrema(imask, "spots/simu_I6_PTmax.txt", "PTmax", patchsize)


  call imap%stack2fig("spots/simu_I6_Tmax.txt", "T", patchsize, output_dir//"stack_simu_deep6_T_onTmax.txt", caption="$T$ on $T_{\max}$", label = "$T (\mu K)$", color_table = "Rainbow")
  call imap%stack2fig("spots/simu_I6_Tmax_QTUTOrient.txt", "T", patchsize, output_dir//"stack_simu_deep6_T_onTmax_Oriented.txt", caption="$T$ on oriented $T_{\max}$", label = "$T (\mu K)$", color_table = "Rainbow")

  call imap%stack2fig("spots/simu_I6_Tmax_QTUTOrient.txt", "Q", patchsize, output_dir//"stack_simu_deep6_Q_onTmax_Oriented.txt", caption="$Q$ on oriented $T_{\max}$", label = "$Q (\mu K)$", color_table = "Rainbow")
  call imap%stack2fig("spots/simu_I6_PTmax.txt", "Q", patchsize, output_dir//"stack_simu_deep6_Q_onPTmax_Oriented.txt" , caption="$Q$ on oriented $P_{T,\max}$", label = "$Q (\mu K)$", color_table = "Rainbow")
  call imap%stack2fig("spots/simu_I6_Tmax.txt", "Qr", patchsize, output_dir//"stack_simu_deep6_Qr_onTmax.txt", caption="$Q_r$ on $T_{\max}$", label = "$Q_r (\mu K)$", color_table = "Rainbow")


  call imap%QU2EB()
  call imap%stack2fig("spots/simu_I6_Tmax.txt", "E", patchsize, output_dir//"stack_simu_deep6_E_onTmax.txt", caption="$E$ on $T_{\max}$", label = "$E (\mu K)$", color_table = "Rainbow")


end program test
