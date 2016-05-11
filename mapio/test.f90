program shells
  use coop_wrapper_firstorder
  use coop_zeta3d_mod
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  implicit none
#include "constants.h"
  type(coop_asy)::fig
  COOP_INT,parameter::lmin = 2, lmax = 600
  COOP_REAL::Cls_Theory(lmin:lmax, 3), Cls_data(lmin:lmax, 3), ells(lmin:lmax)
  type(coop_file)::fp
  COOP_INT::l
  call fp%open_skip_comments("lcdm_scalCls.dat")
  do l = lmin, lmax
     read(fp%unit, *) ells(l), Cls_Theory(l, :)
  enddo
  call fp%close()
  call fp%open_skip_comments("zetaproj/lcdm2_600_Cls.txt")
  do l = lmin, lmax
     read(fp%unit, *) ells(l), Cls_data(l, :)
  enddo
  call fp%close()
  call fig%open("Cl_TT.txt")
  call fig%init(xlabel="$\ell$", ylabel="$\ell(\ell+1)C_l^{TT}/(2\pi)$", xlog = .false.)
  call fig%band(x=ells, ylower=Cls_theory(:,1)*(1.d0-2.d0/sqrt(ells+0.5d0)), yupper =Cls_theory(:,1)*(1.d0+2.d0/sqrt(ells+0.5d0)), colorfill = "RGB:150:150:255", linecolor="invisible")
  call fig%band(x=ells, ylower=Cls_theory(:,1)*(1.d0-1.d0/sqrt(ells+0.5d0)), yupper=Cls_theory(:,1)*(1.d0+1.d0/sqrt(ells+0.5d0)), colorfill = "RGB:100:100:180", linecolor="invisible")
  call fig%dots(ells, Cls_data(:,1), color="RGB:50:50:50")
  call fig%plot(ells, Cls_theory(:,1), linewidth =2., color="red")
  call fig%close()


  call fig%open("Cl_EE.txt")
  call fig%init(xlabel="$\ell$", ylabel="$\ell(\ell+1)C_l^{EE}/(2\pi)$", xlog = .false.)
  call fig%band(x=ells, ylower=Cls_theory(:,2)*(1.d0-2.d0/sqrt(ells+0.5d0)), yupper =Cls_theory(:,2)*(1.d0+2.d0/sqrt(ells+0.5d0)), colorfill = "RGB:150:150:255", linecolor="invisible")
  call fig%band(x=ells, ylower=Cls_theory(:,2)*(1.d0-1.d0/sqrt(ells+0.5d0)), yupper=Cls_theory(:,2)*(1.d0+1.d0/sqrt(ells+0.5d0)), colorfill = "RGB:100:100:180", linecolor="invisible")
  call fig%dots(ells, Cls_data(:,2), color="RGB:50:50:50")
  call fig%plot(ells, Cls_theory(:,2), linewidth =2., color="red")
  call fig%close()


  call fig%open("Cl_TE.txt")
  call fig%init(xlabel="$\ell$", ylabel="$\ell(\ell+1)C_l^{TE}/(2\pi)$", xlog = .false.)
  call fig%band(x=ells, ylower=Cls_theory(:,3)-sqrt((Cls_theory(:,1)*Cls_theory(:,2)+cls_theory(:,3)**2)/(2.d0*ells+1.d0)), yupper =Cls_theory(:,3)+sqrt((Cls_theory(:,1)*Cls_theory(:,2)+cls_theory(:,3)**2)/(2.d0*ells+1.d0)), colorfill = "RGB:150:150:255", linecolor="invisible")
  call fig%band(x=ells, ylower=Cls_theory(:,3)-2.d0*sqrt((Cls_theory(:,1)*Cls_theory(:,2)+cls_theory(:,3)**2)/(2.d0*ells+1.d0)), yupper=Cls_theory(:,3)+2.d0*sqrt((Cls_theory(:,1)*Cls_theory(:,2)+cls_theory(:,3)**2)/(2.d0*ells+1.d0)), colorfill = "RGB:100:100:180", linecolor="invisible")
  call fig%dots(ells, Cls_data(:,3), color="RGB:50:50:50")
  call fig%plot(ells, Cls_theory(:,3), linewidth =2., color="red")
  call fig%close()


end program shells
