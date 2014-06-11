program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"

  character(LEN=*),parameter::mapdir = "act/"
  character(LEN=*),parameter::fitsfile = mapdir//"I50.fits"
  character(LEN=*),parameter::sourcefile = mapdir//"beam50.fits"
  character(LEN=*),parameter::maskfile = mapdir//"W_I50.fits"
  type(coop_fits_image_cea)::cf, mask, src
  type(coop_asy)::asy
  integer, parameter::n=300
  integer ix, iy
  COOP_REAL map(n, n)
  COOP_REAL, parameter::smooth_scale = coop_SI_arcmin * 2.
  call cf%open(fitsfile)
  call cf%get_data()
  call mask%open(maskfile)
  call mask%get_data()
  call src%open(sourcefile)
  call src%get_data
  cf%image = (cf%image)*mask%image/maxval(mask%image)
  print*, maxval(cf%image), minval(cf%image)
  call cf%regularize(0.002d0)
  call cf%filter(lmin = 200, lmax  = 3000)
  print*, maxval(cf%image), minval(cf%image)

  ix = cf%nside(1)/2
  iy = cf%nside(2)/2
  call cf%cut(ix, iy, map)
  call asy%open("cutmap.txt")
  call asy%init(xlabel = "$x$", ylabel  ="$y$", width = 8., height=5.)
  call coop_asy_density(asy, map, 0.d0, dble(n*0.0083*coop_SI_degree), 0.d0, dble(n*0.0083*coop_SI_degree), label = "$I(\mu K)$")
  call asy%close()

!!$  call src%cut(ix, iy, map)
!!$  call asy%open("cutsrc.txt")
!!$  call asy%init(xlabel = "$x$", ylabel  ="$y$", width = 8., height=5.)
!!$  call coop_asy_density(asy, map, 0., real(n*0.0083*coop_SI_degree), 0., real(n*0.0083*coop_SI_degree), label = "$I(\mu K)$")
!!$  call asy%close()

!!$

  call cf%get_flatmap(smooth_scale)
  call asy%open("map.txt")
  call asy%init(xlabel = "$x$", ylabel  ="$y$", width = 8., height=5.)
  call coop_asy_density(asy, cf%smooth_image, cf%xmin, cf%xmax, cf%ymin, cf%ymax, label = "$I(\mu K)$")
  call asy%close()
!!$
!!$  call mask%get_flatmap(smooth_scale)
!!$  call asy%open("mask.txt")
!!$  call asy%init(xlabel = "$x$", ylabel  ="$y$", width = 8., height=5.)
!!$  call coop_asy_density(asy, mask%smooth_image, mask%xmin, mask%xmax, mask%ymin, mask%ymax, label = "weight")
!!$  call asy%close()
!!$
!!$  call src%get_flatmap(smooth_scale)
!!$  call asy%open("src.txt")
!!$  call asy%init(xlabel = "$x$", ylabel  ="$y$", width = 8., height=5.)
!!$  call coop_asy_density(asy, src%smooth_image, src%xmin, src%xmax, src%ymin, src%ymax, label = "source")
!!$  call asy%close()
!!$
end program test
