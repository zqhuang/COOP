program test
#include "constants.h"      
  use coop_wrapper_utils
  implicit none  
  COOP_INT,parameter::n=256
  COOP_REAL::H0(n), P(n)
  COOP_INT::lbd, rbd
  type(coop_asy)::fig
  call coop_set_uniform(n, H0, 63.d0, 78.d0)
  call fig%open("H0posterior.txt")
  call fig%init(width=6., height=4.5, xlabel="$H_0\ (\mathrm{km\,s^{-1}Mpc^{-1}})$", ylabel="$P/P_{\max}$", ymin = 0., ymax=1.05, xmin = 63., xmax = 78.)
  P = exp(-((H0-66.55)/0.89)**2/2.)
  call get_bds()  
  call fig%plot(x= H0, y=P, linewidth=2., color="black", legend="Planck + Wavelet", linetype="solid")

  P = exp(-((H0-67.15)/0.60)**2/2.)
  call get_bds()
  call fig%plot(x= H0, y=P, linewidth=2., color=coop_asy_rgb_color(0.1, 0.2, 0.7), legend="Planck + 12-knots-spline ", linetype="dotted")

  
  P = exp(-((H0-67.36)/0.54)**2/2.)
  call get_bds()  
  call fig%plot(x=H0(lbd:rbd), y=P(lbd:rbd),legend="Planck + $\Lambda$CDM", color=coop_asy_rgb_color(0.2, 0.7, 0.2), linewidth=2., linetype="dashed")


  P = exp(-((H0-73.82)/1.1)**2/2.)
  call get_bds()  
  call fig%plot(x=H0(lbd:rbd), y=P(lbd:rbd),legend="SH0ES + H0LiCOW", color="orange", linewidth=2., linetype="dotdashed")
  
  call fig%top_legend(2)
  call fig%close()
contains
  subroutine get_bds()
    lbd=1
    rbd = n
    do while(P(lbd).lt. 2.d-4)
       lbd = lbd+1
    enddo
    do while(P(rbd).lt. 2.d-4)
       rbd = rbd-1
    enddo
       
  end subroutine get_bds
end program test



