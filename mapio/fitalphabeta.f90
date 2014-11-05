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
  integer, parameter::lmin = 50 
  integer, parameter::lmax = 1500
  type(coop_asy)::fig
  type(coop_file)::fp
  integer l
  integer,parameter::m = 3
  COOP_REAL::scal(m), y(lmin:lmax), tpls(lmin:lmax, m)
  COOP_REAL ells(0:lmax), factor(0:lmax), lnl(0:lmax)
  COOP_REAL,dimension(0:lmax, 1:3) ::Cls_noise, Cls_cmb, Cls_total, Cls_data
  
  call fp%open(coop_inputArgs(1), "r")
  do l=0, lmax
     factor(l) = l*(l+1.d0)/coop_2pi*1.e12
     lnl(l) = log(l+1.d-15)
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
     Cls_total(l, :) = Cls_total(l, :)
  enddo
  call fp%close()

  call fp%open(coop_inputArgs(4), "r")
  do l=0, lmax
     read(fp%unit, *) ells(l), Cls_data(l, :)
  enddo
  call fp%close()

  y(lmin:lmax) = cls_data(lmin:lmax,1)*factor(lmin:lmax)
  tpls(lmin:lmax,1)=cls_cmb(lmin:lmax, 1)*factor(lmin:lmax)
  tpls(lmin:lmax,2)=cls_cmb(lmin:lmax, 1)*factor(lmin:lmax)*ells(lmin:lmax)
  tpls(lmin:lmax,3)=cls_noise(lmin:lmax, 1)*factor(lmin:lmax)
  call coop_fit_template(lmax-lmin+1, m, y, tpls, scal)
  print *, scal

  call fig%open(coop_inputArgs(5))
  call fig%init(xlabel = "$\ell$", ylabel = "$\frac{\ell (\ell + 1) }{2\pi} C_{\ell}^{TE}$")
  call fig%curve(x = ells(lmin:lmax), y = Cls_data(lmin:lmax, 1)*factor(lmin:lmax), legend = "data", color = "red")
  call fig%curve(x = ells(lmin:lmax), y = (Cls_cmb(lmin:lmax, 1)*scal(1)+Cls_noise(lmin:lmax, 1)*(scal(2)+scal(3)*lnl(lmin:lmax)))*factor(lmin:lmax), legend = "FFP8 (cmb + noise) $C_\ell$", color = "blue" , linetype="dotted")
  call fig%legend(0.5, 0.9)
  call fig%close()

  
end program test
