program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  integer,parameter::lmin = 220
  integer,parameter::lmax = 2000
  integer,parameter::irepeat = 1
  COOP_REAL, parameter::reg_limit = 0.003
  character(LEN=*),parameter::mapdir = "act15/"
  character(LEN=*),parameter::Ifile = mapdir//"dataCoadd_I_4.fits"
  character(LEN=*),parameter::Qfile = mapdir//"dataCoadd_Q_4.fits"
  character(LEN=*),parameter::Ufile = mapdir//"dataCoadd_U_4.fits"
  character(LEN=*),parameter::Imaskfile = mapdir//"weightMap_4.fits"
  type(coop_fits_image_cea)::imap, umap, qmap, imask, cutmask
  type(coop_asy)::asy
  integer, parameter::n=300
  integer ix, iy, i, l
  type(coop_file) fp
  COOP_REAL, parameter::patchsize = 90.d0*coop_SI_arcmin
  COOP_UNKNOWN_STRING,parameter::output_dir = "ACTstacking/"

  COOP_REAL map(n, n), Cls(lmin:lmax)
  COOP_REAL, parameter::smooth_scale = coop_SI_arcmin * 5.d0
  call imap%open(Ifile)
  call imap%header%print()
!  call imask%open(Imaskfile)
  call qmap%open(Qfile)
  call umap%open(Ufile)
  print*, maxval(imap%image), minval(imap%image)
  print*, maxval(qmap%image), minval(qmap%image)
  print*, maxval(umap%image), minval(umap%image)  
  call imap%regularize(reg_limit)
  call qmap%regularize(reg_limit)
  call umap%regularize(reg_limit)
  print*, maxval(imap%image), minval(imap%image)
  print*, maxval(qmap%image), minval(qmap%image)
  print*, maxval(umap%image), minval(umap%image)  
  call imap%get_flatmap(smooth_scale)
  imask = imap
  imask%smooth_image = 1.
  
  where(qmap%image**2+umap%image**2 .gt. 2.e5)
     qmap%image = 0.
     umap%image = 0.
  end where

  call qmap%get_flatmap(smooth_scale)
  call umap%get_flatmap(smooth_scale)
  call imap%get_QU(qmap, umap)
  call imap%smooth_flat(lmin = lmin, lmax = lmax)
  call qmap%smooth_flat(lmin = lmin, lmax = lmax)
  call umap%smooth_flat(lmin = lmin, lmax = lmax)  

  call coop_random_init()
  call imap%find_extrema(imask, "spots/act_Tmax.txt", "Tmax", patchsize, irepeat)
  call imap%stack2fig("spots/act_Tmax.txt", "T", patchsize, output_dir//"act_T_onTmax.txt", caption="$T$ on $T_{\max}$", label = "$T (\mu K)$", color_table = "Rainbow")
  call system("../utils/fasy.sh "//output_dir//"act_T_onTmax.txt")  
  call imap%stack2fig("spots/act_Tmax.txt", "Qr", patchsize, output_dir//"act_Qr_onTmax.txt", caption="$Q_r$ on $T_{\max}$", label = "$Q_r (\mu K)$", color_table = "Rainbow")
  call system("../utils/fasy.sh "//output_dir//"act_Qr_onTmax.txt")

  call imap%stack2fig("spots/act_Tmax.txt", "Q", patchsize, output_dir//"act_Q_onTmax.txt", caption="$Q$ on $T_{\max}$", label = "$Q (\mu K)$", color_table = "Rainbow")
  call system("../utils/fasy.sh "//output_dir//"act_Q_onTmax.txt")
  

!!simulations

!!$  call imap%simulate_flat(lmin = lmin,lmax = lmax, cls_file = "cls.dat")
!!$  call imap%find_extrema(cutmask, "spots/simu_I"//mapid//"_Tmax.txt", "Tmax", patchsize, irepeat)
!!$
!!$  call imap%stack2fig("spots/simu_I"//mapid//"_Tmax.txt", "Qr", patchsize, output_dir//"simu_Qr_onTmax.txt", caption="$Q_r$ on $T_{\max}$", label = "$Q_r (\mu K)$", color_table = "Rainbow")



end program test
