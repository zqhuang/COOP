program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  implicit none
#include "constants.h"
  integer, parameter::lmax = 2048
  type(coop_asy)::fig
  type(coop_file)::fp
  integer l
  COOP_REAL ells(0:lmax), factor(0:lmax)
  COOP_REAL,dimension(0:lmax, 1:3) ::Cls_noise, Cls_cmb, Cls_total, Cls_data
  call fp%open(coop_inputArgs(1), "r")
  do l=0, lmax
     factor(l) = l*(l+1.d0)/coop_2pi*1.e12
     read(fp%unit, *) ells(l), Cls_cmb(l, :)
  enddo
  call fp%close()
  
  call fp%open(coop_inputArgs(2), "r")
  do l=0, lmax
     read(fp%unit, *) ells(l), Cls_noise(l, :)
  enddo
  call fp%close()

  call fp%open(coop_inputArgs(3), "r")
  do l=0, lmax
     read(fp%unit, *) ells(l), Cls_total(l, :)
  enddo
  call fp%close()

  call fp%open(coop_inputArgs(4), "r")
  do l=0, lmax
     read(fp%unit, *) ells(l), Cls_data(l, :)
  enddo
  call fp%close()

  call fig%open(coop_inputArgs(5))
  call fig%init(xlabel = "$\ell$", ylabel = "$\frac{\ell (\ell + 1) }{2\pi} C_{\ell}^{TE}$")
  call fig%curve(x = ells(2:lmax), y = Cls_data(2:lmax, 1)*factor(2:lmax), legend = "data", color = "red")
  call fig%curve(x = ells(2:lmax), y = Cls_total(2:lmax, 1)*factor(2:lmax), legend = "FFP8 (cmb + noise) $C_\ell$", color = "blue" , linetype="dotted")
  call fig%curve(x = ells(2:lmax), y = Cls_noise(2:lmax, 1)*factor(2:lmax), color = "gray", linewidth = 0.5, legend = "FFP8 noise")
  call fig%curve(x = ells(2:lmax), y = Cls_cmb(2:lmax, 1)*factor(2:lmax), color = "cyan", linewidth = 0.5, legend = "FFP8 cmb")
  call fig%curve(x = ells(2:lmax), y = (Cls_cmb(2:lmax, 1)+Cls_noise(2:lmax, 1))*factor(2:lmax), legend = "FFP8 noise $C_\ell$ + cmb $C_\ell$", linetype="dashed")
  call fig%legend(0.5, 0.9)
  call fig%close()

  
end program test
