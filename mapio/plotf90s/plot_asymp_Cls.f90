program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
#ifdef HAS_HEALPIX
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  use udgrade_nr
  use coord_v_convert,only:coordsys2euler_zyz
#endif  
  implicit none
#include "constants.h"
  COOP_INT, parameter::lmax = 30
  type(coop_file)::fp
  type(coop_asy)::fig
  COOP_SINGLE, parameter::linewidth = 0.5
  COOP_STRING::line
  COOP_INT::l, i , ic
  COOP_REAL::tmp
  COOP_REAL::ells(2:lmax), Cls(2:lmax), Cls_fiducial(2:lmax), delta_Cls(2:lmax)
  call coop_set_uniform(lmax-1, ells, 2.d0, dble(lmax))
  call fp%open("clsout/clsout200.dat","r")
  do l=2, lmax
     read(fp%unit, *) i, Cls_fiducial(l), Cls(l)
     delta_Cls(l) = Cls_fiducial(l)/sqrt(2.d0*l+1.d0)
  enddo
  call fig%open("asymptotic_Cl.txt")
  call fig%init(xlabel = "$\ell$", ylabel = "$\frac{\ell(\ell +1)}{2\pi}C_l \, (\mu K^2)$", width = 4.8, height = 4., xmin = 1.5, xmax = 30.5, ymin = 0., ymax = 2000., doclip=.true.)

  call fig%band(x = ells, ylower = Cls_fiducial - delta_Cls, yupper = Cls_fiducial+delta_Cls, colorfill = coop_asy_gray_color(0.6), smooth = .false., linecolor = "invisible")
  call fig%curve(x = ells, y = Cls_fiducial, color = "black", linetype = "solid", linewidth = 2., legend = "fiducial")  
  call fig%curve(x = ells, y = Cls, color = "red", linetype = "dotted", linewidth = 2., legend = "reconstructed")  
  
  call fig%legend(0.1, 0.92, 2, box = .false.)
  
end program test
