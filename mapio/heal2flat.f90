program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_file)::fp
  type(coop_healpix_maps)::hpT, hpQ, hpU, hpQr, hpUr
  COOP_INT i, pix
  COOP_REAL f(5)
  type(coop_healpix_disc) disc
  type(coop_healpix_patch) patch
  COOP_INT,parameter::n = 36
  COOP_REAL,parameter::r_degree  = 2.d0
  COOP_REAL,parameter::dr = 2.d0*sin(r_degree*coop_SI_degree/2.d0)/n
  call hpT%init(nside = 1024, nmaps = 1, spin = (/ 0 /) )
  hpQ = hpT
  hpU = hpT
  hpQr = hpT
  hpUr = hpT
  call fp%open("raul/patch_commander.dat", "r")
  do i = 0, 23979
     read(fp%unit, *) pix, f
     f = f *1.d6
     if(pix .ge. 0 .and. pix .le. hpT%npix - 1)then
        hpT%map(pix, 1) = f(1)
        hpQ%map(pix, 1) = f(2)
        hpU%map(pix, 1) = f(3)
        hpQr%map(pix, 1) = f(4)
        hpUr%map(pix, 1) = f(5)
     else
        print*, pix
        stop
     endif
  enddo
  call fp%close()
  coop_healpix_patch_default_want_label  = .true.
  coop_healpix_patch_default_figure_width = 3.5

  call hpT%get_disc(0, disc)

  call patch%init("T", n, dr)
  call hpT%fetch_patch(disc, coop_pio2, patch)
  call patch%plot(1, "raul/raulT.txt")
  call hpQ%fetch_patch(disc, coop_pio2, patch)
  call patch%plot(1, "raul/raulQ.txt")
  call hpU%fetch_patch(disc, coop_pio2, patch)
  call patch%plot(1, "raul/raulU.txt")
  call hpQr%fetch_patch(disc, coop_pio2, patch)
  call patch%plot(1, "raul/raulQr.txt")
  call hpUr%fetch_patch(disc, coop_pio2, patch)
  call patch%plot(1, "raul/raulUr.txt")
  
end program test
