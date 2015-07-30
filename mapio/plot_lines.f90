program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  implicit none
#include "constants.h"
  COOP_INT, parameter::n = 21, im = 1
  type(coop_asy)::fig
  COOP_REAL::nu(n), f(n)
  call fig%open("minkowski"//COOP_STR_OF(im)//".txt")
  select case(im)
  case(0)
     call fig%init(xlabel = "$\nu$", ylabel = "$dA/d\nu \ln[\sqrt{2\pi}e^{\nu^2/2}dA/d\nu]$", ymax = 0.1, ymin = -0.1)
  case(1)
     call fig%init(xlabel = "$\nu$", ylabel = "$ \sqrt{\frac{4}{\pi}} \frac{C \sigma_0} {\sigma_1} \ln[ 2\sqrt{2} e^{\nu^2}  \frac{C \sigma_0} {\sigma_1} ] $", ymax = 0.05, ymin = -0.1)     
  end select
  call load_data("dust_15a_V"//COOP_STR_OF(im)//".txt")
  call fig%curve(nu, f, color = "red", linetype = "solid", linewidth = 1.5, legend = "dust 15'")
  call load_data("gauss_15a_V"//COOP_STR_OF(im)//".txt")
  call fig%curve(nu, f, color = "red", linetype = "dotted", linewidth = 1.5, legend = "Gaussian sim. 15'")
  call load_data("dust_30a_V"//COOP_STR_OF(im)//".txt")
  call fig%curve(nu, f, color = "blue", linetype = "solid", linewidth = 1.5, legend = "dust 30a'")
  call load_data("gauss_30a_V"//COOP_STR_OF(im)//".txt")
  call fig%curve(nu, f, color = "blue", linetype = "dotted", linewidth = 1.5, legend = "Gaussian sim. 30'")
  call load_data("dust_60a_V"//COOP_STR_OF(im)//".txt")
  call fig%curve(nu, f, color = "black", linetype = "solid", linewidth = 1.5, legend = "dust 60a'")
  call load_data("gauss_60a_V"//COOP_STR_OF(im)//".txt")
  call fig%curve(nu, f, color = "black", linetype = "dotted", linewidth = 1.5, legend = "Gaussian sim. 60'")
  f = 0.
  call fig%curve(nu, f, color="gray", linetype = "dashed", linewidth = 0.5)
  call fig%legend_advance(0.1, 0.25,  box_color = "invisible", xmargin = 0., ymargin = 0., linelength = 0.8, hskip=0.7, vskip = 0.85, cols = 2)
  call system("../utils/fasy.sh minkowski"//COOP_STR_OF(im)//".txt")
  call fig%close()

contains
  subroutine load_data(fname)
    type(coop_file)::fp
    COOP_UNKNOWN_STRING::fname
    COOP_INT::i
    call fp%open(fname, "r")
    do i=1, n
       read(fp%unit, *) nu(i), f(i)
    enddo
    call fp%close()
  end subroutine load_data
end program test
