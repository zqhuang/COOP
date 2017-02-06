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
  COOP_UNKNOWN_STRING,parameter::mode = "spn"
  COOP_UNKNOWN_STRING,parameter::unit = "muK"
  COOP_INT,parameter::npt = 41, nsim  = 30
  type(coop_file)::fdata, fsim
  COOP_STRING::legend, simlegend
  COOP_INT::i, ipt
  COOP_REAL::ydata(npt), ysim(npt), ys(npt, nsim), x(npt), xtmp,  ysim2(npt)
  COOP_SINGLE::ypos, xpos, ypos2
  type(coop_asy)::fig
  xpos = 0.3
  select case(trim(field))
  case("T")
     ypos = 0.7
     ypos2 = 0.5
  case("E")
     ypos = 0.5
     ypos2 = 0.3
  case("B")
     ypos = 0.9
     ypos2 = 0.3
  case default
     print*, "field = "//trim(field)
     stop "Unknown field"
  end select
  select case(trim(mode))
  case("spn")
     call fdata%open("actpics/nF0_cutx90y50_act_on_cutx90y50_act_l"//COOP_STR_OF(ell)//"_5a_randrot_"//trim(field)//"onThotcoldnu0_m0.dat")
     legend = "ACT $"//trim(field)//"$ on ACT $T$"
     simlegend = "$\Lambda$CDM + noise sim. $1\sigma$, $2\sigma$"
  case("s")
     call fdata%open("actpics/nF0_cutx90y50_act_on_cutx90y50_planck_l"//COOP_STR_OF(ell)//"_5a_randrot_"//trim(field)//"onThotcoldnu0_m0.dat")
     legend = "ACT $"//trim(field)//"$ on PLANCK $T$"     
     simlegend = "$\Lambda$CDM"
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
  call fig%init(xlabel="$r$ [deg]", ylabel = "$"//trim(field)//"(\mu K)$")

  do i = 1, nsim
     call fsim%open("actpics/nF0_cutx90y50_"//trim(mode)//COOP_STR_OF(i)//"_on_cutx90y50_"//trim(mode)//COOP_STR_OF(i)//"_l"//COOP_STR_OF(ell)//"_5a_randrot_"//trim(field)//"onThotcoldnu0_m0.dat")
     do ipt = 1, npt
        read(fsim%unit, *) xtmp, ys(ipt, i)
        if(abs(xtmp - x(ipt)).gt. 1.d-4) stop "xaxes do not have the same resolution"
     enddo
     call fsim%close()
     ysim  = ysim + ys(:, i)
     ysim2 = ysim2 + ys(:, i)**2
  enddo
  ysim = ysim/nsim
  ysim2 = sqrt(max(ysim2/nsim - ysim**2, 0.d0))
  x = x/coop_SI_degree
  call fig%band(x=x, ylower = ysim - ysim2*2, yupper = ysim+ysim2*2, colorfill = coop_asy_gray_color(0.75), smooth = .false., linecolor = coop_asy_gray_color(0.75) , linewidth = 1.)
  
  call fig%band(x=x, ylower = ysim - ysim2, yupper = ysim+ysim2, colorfill =  coop_asy_gray_color(0.6), smooth = .false., linecolor =  coop_asy_gray_color(0.6) , linewidth = 1., legend=trim(simlegend))
!!$  do i = 1, nsim, 4
!!$     if(i .eq. 1)then
!!$        call fig%plot(x, ys(:, i), color="blue", linetype="dashed", linewidth=1., legend="simulations")
!!$     else
!!$        call fig%plot(x, ys(:, i), color="blue", linetype="dashed", linewidth=1.)
!!$     endif
!!$  enddo

  call fig%plot(x , ydata, color = "red", linetype = "dotted", linewidth=2.5, legend=trim(legend))
  call fig%legend(xpos, ypos)
  call fig%label("$\ell_{\min} = "//COOP_STR_OF(ell)//"$, FWHM=5 arcmin", xpos, ypos2)
  call fig%label("$\ell_{x,\min} = 90,\ \ell_{y,\min}=50$", xpos, ypos2-0.1)
  call fig%close()
end program shells
