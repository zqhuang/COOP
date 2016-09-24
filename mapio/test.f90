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
  COOP_INT,parameter::ell = 350
  COOP_UNKNOWN_STRING,parameter::field = "T"
  COOP_UNKNOWN_STRING,parameter::mode = "s"
  COOP_UNKNOWN_STRING,parameter::unit = "muK"
  COOP_INT,parameter::npt = 41, nsim  =30
  type(coop_file)::fdata, fsim
  COOP_STRING::legend
  COOP_INT::i, ipt
  COOP_REAL::ydata(npt), ysim(npt), x(npt), xtmp, ytmp(npt), ysim2(npt)
  type(coop_asy)::fig
  select case(mode)
  case("spn")
     call fdata%open("actpics/nF0_cutx90y50_act_on_cutx90y50_act_l"//COOP_STR_OF(ell)//"_5a_randrot_"//trim(field)//"onThotcoldnu0_m0.dat")
     legend = "ACT $"//trim(field)//"$ on ACT $T$"
  case("s")
     call fdata%open("actpics/nF0_cutx90y50_act_on_cutx90y50_planck_l"//COOP_STR_OF(ell)//"_5a_randrot_"//trim(field)//"onThotcoldnu0_m0.dat")
     legend = "ACT $"//trim(field)//"$ on PLANCK $T$"
  case default
     stop "unknown mode: the mode must be s or spn"
  end select
  do ipt = 1, npt
     read(fdata%unit,*) x(ipt), ydata(ipt)
  enddo
  call fdata%close()
  ysim = 0.d0
  ysim2 = 0.d0
  call fig%open("l"//COOP_STR_OF(ell)//"_"//trim(mode)//"_"//trim(field)//".txt")
  call fig%init(xlabel="$\varpi$", ylabel = "$"//trim(field)//"(\mu K)$")

  do i = 1, nsim
     call fsim%open("actpics/nF0_cutx90y50_"//trim(mode)//COOP_STR_OF(i)//"_on_cutx90y50_"//trim(mode)//COOP_STR_OF(i)//"_l"//COOP_STR_OF(ell)//"_5a_randrot_"//trim(field)//"onThotcoldnu0_m0.dat")
     do ipt = 1, npt
        read(fsim%unit, *) xtmp, ytmp(ipt)
        if(abs(xtmp - x(ipt)).gt. 1.d-4) stop "xaxes do not have the same resolution"
     enddo
     call fsim%close()
     ysim  = ysim + ytmp
     ysim2 = ysim2 + ytmp**2
  enddo
  ysim = ysim/nsim
  ysim2 = sqrt(max(ysim2/nsim - ysim**2, 0.d0))
  call fig%band(x=x, ylower = ysim - ysim2*2, yupper = ysim+ysim2*2, colorfill = coop_asy_gray_color(0.7), smooth = .false., linecolor = coop_asy_gray_color(0.7) , linewidth = 1.)
  call fig%band(x=x, ylower = ysim - ysim2, yupper = ysim+ysim2, colorfill = "gray", smooth = .false., linecolor = "gray" , linewidth = 1.)
  call fig%plot(x , ydata, color = "red", linetype = "solid", linewidth=1.5, legend=trim(legend))
  call fig%legend(0.3, 0.6)
  call fig%label("$\ell_{\min} = "//COOP_STR_OF(ell)//"$, FWHM=5 arcmin", 0.3, 0.5)
  call fig%label("$\ell_{x,\min} = 90,\ \ell_{y,\min}=50$", 0.3, 0.4)
  call fig%close()
end program shells
