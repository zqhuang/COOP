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
  COOP_REAL::tmp, tmp2
  COOP_REAL::ells(2:lmax), Cls(2:lmax), Cls_fiducial(2:lmax), delta_Cls(2:lmax), delta_mCls(2:lmax)
  call coop_set_uniform(lmax-1, ells, 2.d0, dble(lmax))
  call fig%open("lowl_Cl.txt")
  call fig%init(xlabel = "$\ell$", ylabel = "$\frac{\ell(\ell +1)}{2\pi}C_l \, (\mu K^2)$", width = 6., height = 5., xmin = 1.5, xmax = 30.5, ymin = 0., ymax = 2200., doclip=.true.)

  call fp%open("plikv17_cls.dat", "r")
  do l = 2, lmax
     read(fp%unit, *) i, tmp, tmp2, Cls(l), delta_Cls(l), delta_mCls(l)
     if(i.ne. l) stop "error"
  enddo
  call fp%close()


  call fig%curve(x = ells, y = Cls, legend="Commander NSIDE = 16; From Likelihood Group", color ="gray", linetype = "solid", linewidth = 2.)  
  call fig%errorbars(x = ells, y = Cls, dy = delta_Cls, dy_minus = delta_Cls, color = "gray_solid_1",  barsize = 5., center_color = "gray_solid_4")
  

  call fp%open("clsout/clsout_NONE.dat", "r")
  read(fp%unit, "(A)")line
  do l = 2, lmax
     read(fp%unit, *) i, Cls_fiducial(l), Cls(l), delta_Cls(l)
     if(i.ne. l) stop "error"
  enddo
  call fp%close()

  call fig%curve(x = ells, y = Cls_fiducial, legend="$\Lambda$CDM", color = "red", linetype = "solid", linewidth = 1.)
  
  call fig%curve(x = ells, y = Cls, legend="Commander NSIDE=256; From Zhiqi; Ignoring noise", color ="blue", linetype = "dotted", linewidth = 2.)  
  

  
  
  call fig%legend(0.1, 0.92, 1, box = .false.)
  
end program test
