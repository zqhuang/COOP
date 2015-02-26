program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  implicit none
#include "constants.h"
  integer,parameter::lmin = 200
  integer,parameter::lmax = 2000
  integer,parameter::irepeat = 1
  COOP_REAL, parameter::reg_limit = 0.01
  character(LEN=*),parameter::mapdir = "act15/"
  character(LEN=*),parameter::Ifile = mapdir//"dataCoadd_I_4.fits"
  character(LEN=*),parameter::Qfile = mapdir//"dataCoadd_Q_4.fits"
  character(LEN=*),parameter::Ufile = mapdir//"dataCoadd_U_4.fits"
  character(LEN=*),parameter::Imaskfile = mapdir//"weightMap_4.fits"
  type(coop_fits_image_cea)::imap, umap, qmap
  type(coop_healpix_maps)::hp, imask, qmask, umask
  call imap%open(Ifile)
  call qmap%open(Qfile)
  call umap%open(Ufile)
  print*, maxval(imap%image), minval(imap%image)
  print*, maxval(qmap%image), minval(qmap%image)
  print*, maxval(umap%image), minval(umap%image)  
  call hp%init(nside = 2048, nmaps = 3, spin = (/ 0, 2, 2 /) )
  call imask%init(nside = 2048, nmaps = 1, spin = (/ 0 /) )
  call qmask%init(nside = 2048, nmaps = 1, spin = (/ 0 /) )
  call umask%init(nside = 2048, nmaps = 1, spin = (/ 0 /) )    
  call imap%convert2healpix(hp, 1, imask, lower = -400., upper = 400.)
  call qmap%convert2healpix(hp, 2, qmask, lower = -200., upper = 200.)
  call umap%convert2healpix(hp, 3, umask, lower = -200., upper = 200.)
  print*, maxval(hp%map(:,1)), minval(hp%map(:,1))
  print*, maxval(hp%map(:,2)), minval(hp%map(:,2))
  print*, maxval(hp%map(:,3)), minval(hp%map(:,3))  
  hp%map(:,2:3) = - hp%map(:, 2:3) !!convert to healpix convention
  qmask%map(:,1) = qmask%map(:,1)*umask%map(:,1)
!  call hp%rotate_coor('C','G')
!  call imask%rotate_coor('C','G')
!  call qmask%rotate_coor('C','G')
  call hp%write("act15/act15_i.fits", index_list = (/ 1 /) )
  call hp%write("act15/act15_pol.fits", index_list = (/ 2, 3 /) )  
  call imask%write("act15/act15_imask.fits")
  call qmask%write("act15/act15_polmask.fits")  
end program test
