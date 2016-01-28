program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  COOP_UNKNOWN_STRING,parameter::mapdir = "act16/"
  COOP_UNKNOWN_STRING,parameter::Ifile = mapdir//"deep56_coadd_I_smooth.fits"
  COOP_UNKNOWN_STRING,parameter::Qfile = mapdir//"deep56_coadd_Q_smooth.fits"
  COOP_UNKNOWN_STRING,parameter::Ufile = mapdir//"deep56_coadd_U_smooth.fits"
  COOP_UNKNOWN_STRING,parameter::I_mask_file = mapdir//"deep56_weight_I_smooth.fits"
  COOP_UNKNOWN_STRING,parameter::Q_mask_file = mapdir//"deep56_weight_Q_smooth.fits"
  COOP_UNKNOWN_STRING,parameter::U_mask_file = mapdir//"deep56_weight_U_smooth.fits"
  type(coop_fits_image_cea)::imap, umap, qmap, Imask, Qmask, UMask
  COOP_REAL, parameter::patchsize = 90.d0*coop_SI_arcmin, smooth_scale=1.5d0*coop_SI_arcmin
  COOP_UNKNOWN_STRING,parameter::output_dir = "ACTstacking/"
  call imap%open(IFile)
  call qmap%open(QFile)
  call Umap%open(UFile)
  call imask%open(i_mask_file)
  call qmask%open(q_mask_file)
  call umask%open(u_mask_file)
  where(imask%image .lt. 10.)
     imask%image = 0.
     imap%image = 0.
  elsewhere
     imask%image= 1.
  end where
  where(qmask%image .lt. 10.)
     qmask%image = 0.
     qmap%image = 0.
  elsewhere
     qmask%image= 1.
  end where
  where(umask%image .lt. 10.)
     umask%image = 0.
     umap%image = 0.
  elsewhere
     umask%image= 1.
  end where
  call imap%get_flatmap(smooth_scale)
  call qmap%get_flatmap(smooth_scale)
  call umap%get_flatmap(smooth_scale)
  call imask%get_flatmap(smooth_scale)
  call qmask%get_flatmap(smooth_scale)
  call umask%get_flatmap(smooth_scale)
  call imap%find_extrema(imask, "spots/act_Hot.txt", "Hot", nu=1.d0)
  call imap%stack2fig("spots/act_Hot.txt", "T", patchsize, output_dir//"act_T_onHotnu1.txt", caption="$T$ on $T_{\max}$", label = "$T (\mu K)$", color_table = "Rainbow")
  call imap%get_QU(qmap,umap)
  call imap%stack2fig("spots/act_Hot.txt", "Qr", patchsize, output_dir//"act_Qr_onHotnu1.txt", caption="$Q_r$ on $T_{\max}$", label = "$Q_r (\mu K)$", color_table = "Rainbow")
  call imap%stack2fig("spots/act_Hot.txt", "Q", patchsize, output_dir//"act_Q_onHotnu1.txt", caption="$Q$ on $T_{\max}$", label = "$Q (\mu K)$", color_table = "Rainbow")

   end program test
