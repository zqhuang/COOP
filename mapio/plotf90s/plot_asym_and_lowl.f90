program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_asy)::fig
  type(coop_file)::fp, fp2
  COOP_INT,parameter::lmax = 30
  COOP_REAL::ells(2:lmax), Dls_theory(2:lmax), Dls_data(2:lmax), Dls_low(2:lmax), Dls_up(2:lmax), Delta_Dls(2:lmax)
  COOP_INT::l
  COOP_STRING::line
  call fp%open("clsout/clsout580.dat", "r")
  do l = 2, lmax
     read(fp%unit, *) ells(l), Dls_theory(l), Dls_data(l)
  enddo
  Dls_low = Dls_theory *(1.- sqrt(1./(ells+0.5)))
  Dls_up = Dls_theory *(1.+ sqrt(1./(ells+0.5)))  
  
  call fp%close()
  call fig%open("asymptotic_Cl.txt")
  call fig%init(xlabel="$\ell$", ylabel = "$\frac{\ell(\ell+1)}{2\pi}C_l (\mu K^2)$", xmin = 0., xmax = 32., ymin = 0., ymax = 2000.)
  call fig%band(ells, Dls_low, Dls_up, trim(coop_asy_gray_color(0.6)), .false., "invisible", "solid", 1.)  
  call fig%curve(ells, Dls_theory, legend = "fiducial", color="black", linewidth=2.)
  call fig%curve(ells, Dls_data, legend = "reconstructed", color="red", linewidth=2., linetype = "dotted")
  call fig%legend(0.4, 0.92)
  call fig%close()

  call fp%open("clsout/clsout_NONE.dat", "r")
  read(fp%unit, "(A)") line
  do l = 2, lmax
     read(fp%unit, *) ells(l), Dls_theory(l), Dls_data(l), delta_Dls(l)
  enddo
  call fp%close()
  Dls_low = Dls_data - delta_Dls
  Dls_up = Dls_data + delta_Dls  
  call fig%open("lowl_Cl.txt")
  call fig%init(xlabel="$\ell$", ylabel = "$\frac{\ell(\ell+1)}{2\pi}C_l (\mu K^2$", xmin = 0., xmax = 32.)
  call fig%band(ells, Dls_low, Dls_up, trim(coop_asy_gray_color(0.6)), .false., "invisible", "solid", 1.)
  call fig%curve(ells, Dls_theory, legend = "$\Lambda$CDM", color="black", linewidth=2.)
  call fig%curve(ells, Dls_data, legend = "commander mask", color="red", linewidth=2., linetype = "dotted")

  call fp%open("clsout/clsout_COLDSPOT.dat", "r")
  read(fp%unit, "(A)") line
  do l = 2, lmax
     read(fp%unit, *) ells(l), Dls_theory(l), Dls_data(l), delta_Dls(l)
  enddo
  call fp%close()

  call fig%curve(ells, Dls_data, legend = "commander mask + cold spot", color="blue", linewidth=2., linetype = "dashed")
  call fig%legend(0.24, 0.23)
  call fig%close()
end program test
