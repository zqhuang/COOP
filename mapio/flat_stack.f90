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
  character(LEN=*),parameter::Qmaskfile = mapdir//"WQ6.fits"
  character(LEN=*),parameter::Umaskfile = mapdir//"WU6.fits"
  type(coop_fits_image_cea)::imap, umap, qmap, imask, qmask, umask
  type(coop_asy)::asy
  integer, parameter::n=300
  integer ix, iy, i, l
  type(coop_file) fp
  COOP_REAL map(n, n), Cls(lmin:lmax)
  COOP_REAL, parameter::smooth_scale = coop_SI_arcmin * 1.5
  call imap%open(Ifile)
  call imask%open(Imaskfile)
  imap%image = (imap%image)*imask%image/maxval(imask%image)
  call imap%regularize(0.002d0)
  call imap%get_flatmap(smooth_scale)
  call imap%smooth_flat(lmin = lmin, lmax = lmax)

  call qmap%open(Qfile)
  call qmask%open(Qmaskfile)
  qmap%image = (qmap%image)*qmask%image/maxval(qmask%image)
  call qmap%regularize(0.002d0)
  call qmap%get_flatmap(smooth_scale)
  call qmap%smooth_flat(lmin = lmin, lmax = lmax)

  call umap%open(Ufile)
  call umask%open(Umaskfile)
  umap%image = (umap%image)*umask%image/maxval(umask%image)
  call umap%regularize(0.002d0)
  call umap%get_flatmap(smooth_scale)
  call umap%smooth_flat(lmin = lmin, lmax = lmax)

  call imap%get_QU(qmap, umap)

  call imap%stack2fig("spots/I6_PTmax.txt", "Q", 30.d0*coop_SI_arcmin, "stackedQr.txt")
  
end program test
