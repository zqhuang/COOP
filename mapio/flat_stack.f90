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



  call imap%stack2fig("spots/I6_Tmax_QTUTOrient.txt", "T", 20.d0*coop_SI_arcmin, "stackedT.txt")
  call imap%stack2fig("spots/I6_Tmax_QTUTOrient.txt", "Q", 20.d0*coop_SI_arcmin, "stackedQ.txt")
  call imap%stack2fig("spots/I6_Tmax.txt", "Qr", 30.d0*coop_SI_arcmin, "stackedQr.txt")

  where(imask%image .lt. 5000)
     imask%image = 0.
  elsewhere
     imask%image = 1.
  end where
  call imask%get_flatmap(smooth_scale)

  call imap%simulate_flat(lmin = lmin,lmax = lmax, cls_file = "cls.txt")
  call imap%find_extrema(imask, "spots/simu_I6_Tmax.txt", "Tmax", 30.d0*coop_SI_arcmin)

  qmap = imap
  call qmap%get_QTUT()
  call qmap%find_extrema(imask, "spots/simu_I6_Tmax_QTUTOrient.txt", "Tmax_QTUTOrient", 30.d0*coop_SI_arcmin)


  call imap%stack2fig("spots/simu_I6_Tmax_QTUTOrient.txt", "T", 20.d0*coop_SI_arcmin, "simu_stackedT.txt")

  call imap%stack2fig("spots/simu_I6_Tmax_QTUTOrient.txt", "Q", 20.d0*coop_SI_arcmin, "simu_stackedQ.txt")
  call imap%stack2fig("spots/simu_I6_Tmax.txt", "Qr", 30.d0*coop_SI_arcmin, "simu_stackedQr.txt")


  call imap%QU2EB()
  call imap%stack2fig("spots/simu_I6_Tmax.txt", "E", 30.d0*coop_SI_arcmin, "simu_stackedE.txt")


  call imap%stack2fig("spots/simu_I6_Tmax.txt", "B", 30.d0*coop_SI_arcmin, "simu_stackedB.txt")
  
end program test
