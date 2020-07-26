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
#endif  
  implicit none
#include "constants.h"
  COOP_INT::i, j
  COOP_INT,parameter::n = 30
  COOP_REAL::r(0:n), tn_cold(0:n), ts_cold(0:n), tn_hot(0:n), ts_hot(0:n)
  type(coop_file)::fp
  type(coop_asy)::fig_cold, fig_hot, fig
  COOP_SINGLE::ymin = -8.
  COOP_SINGLE::ymax = 10.
  call fig_cold%open("pwasym_stack_cold.txt")
  call fig_hot%open("pwasym_stack_hot.txt")
  call fig%open("pwasym_stack_hotcold.txt")  
  call fig_cold%init(xlabel = "$\varpi$", ylabel = "$(T_N - T_S) [\mu K]$", ymin = ymin, ymax = ymax, caption = "on $\nu \ge 0.5$ cold pixels; FWHM $2^{\circ}$; high-pass $\ell_{\min} = 10$")
  call fig_hot%init(xlabel = "$\varpi$", ylabel = "$(T_N - T_S) [\mu K]$",  ymin = ymin, ymax = ymax,  caption = "on $\nu\ge 0.5$ hot pixels; FWHM $2^{\circ}$; high-pass $\ell_{\min} = 10$")
  call fig%init(xlabel = "$\varpi$", ylabel = "$(T_N - T_S) [\mu K]$",  ymin = ymin, ymax = ymax, caption = " hot - cold; $\nu\ge 0.5$ pixels; FWHM $2^{\circ}$; high-pass $\ell_{\min} = 10$")    
  
  call fp%open("stacked/nasym_120a_cold_nu0p5_m0.dat", "r")
  do i=0, n
     read(fp%unit, *) r(i), tn_cold(i)
  enddo
  call fp%close()
  call fp%open("stacked/sasym_120a_cold_nu0p5_m0.dat", "r")
  do i=0, n
     read(fp%unit, *) r(i), ts_cold(i)
  enddo
  call fp%close()
  call fp%open("stacked/nasym_120a_hot_nu0p5_m0.dat", "r")
  do i=0, n
     read(fp%unit, *) r(i), tn_hot(i)
  enddo
  call fp%close()
  call fp%open("stacked/sasym_120a_hot_nu0p5_m0.dat", "r")
  do i=0, n
     read(fp%unit, *) r(i), ts_hot(i)
  enddo
  call fp%close()
  
  call fig_hot%curve(x = r, y = (tn_hot - ts_hot), legend = "Planck", color = "red", linewidth = 2.)
  call fig_cold%curve(x = r, y = (tn_cold - ts_cold), legend = "Planck", color = "red", linewidth = 2.)  
  call fig%curve(x = r, y = (tn_hot-tn_cold - ts_hot+ts_cold), legend = "Planck", color = "red", linewidth = 2.)    
  
  do j = 0, 99
     call fp%open("stacked/ffp8_"//COOP_STR_OF(j)//"_nasym_120a_cold_nu0p5_m0.dat", "r")
     do i=0, n
        read(fp%unit, *) r(i), tn_cold(i)
     enddo
     call fp%close()
     call fp%open("stacked/ffp8_"//COOP_STR_OF(j)//"_sasym_120a_cold_nu0p5_m0.dat", "r")
     do i=0, n
        read(fp%unit, *) r(i), ts_cold(i)
     enddo
     call fp%close()
     call fp%open("stacked/ffp8_"//COOP_STR_OF(j)//"_nasym_120a_hot_nu0p5_m0.dat", "r")
     do i=0, n
        read(fp%unit, *) r(i), tn_hot(i)
     enddo
     call fp%close()
     call fp%open("stacked/ffp8_"//COOP_STR_OF(j)//"_sasym_120a_hot_nu0p5_m0.dat", "r")
     do i=0, n
        read(fp%unit, *) r(i), ts_hot(i)
     enddo
     call fp%close()
     if(j.eq.0)then
        call fig_hot%curve(x = r, y = 2.d0*(tn_hot - ts_hot)/(ts_hot+tn_hot), legend = "FFP8", color = "gray", linewidth = 1.2, linetype = "dotted")
        call fig_cold%curve(x = r, y = 2.d0*(tn_cold - ts_cold)/(ts_cold+tn_cold), legend = "FFP8", color = "gray", linewidth = 1.2, linetype = "dotted")  
        call fig%curve(x = r, y = 2.d0*(tn_hot-tn_cold - ts_hot+ts_cold)/(tn_hot-tn_cold+ts_hot - ts_cold), legend = "FFP8", color = "gray", linewidth = 1.2, linetype = "dotted")    
     else
        call fig_hot%curve(x = r, y = (tn_hot - ts_hot), color = "gray", linewidth = 1.2, linetype = "dotted")
        call fig_cold%curve(x = r, y = (tn_cold - ts_cold),  color = "gray", linewidth = 1.2, linetype = "dotted")  
        call fig%curve(x = r, y = (tn_hot-tn_cold - ts_hot+ts_cold), color = "gray", linewidth = 1.2, linetype = "dotted")    
     endif
  enddo
  call fig%legend(0.1, 0.2)
  call fig_hot%legend(0.1, 0.9)
  call fig_cold%legend(0.1, 0.9)  
  call fig%close()
  call fig_hot%close()
  call fig_cold%close()
  call system("../utils/fasy.sh pwasym_stack_hotcold.txt")
  call system("../utils/fasy.sh pwasym_stack_hot.txt")
  call system("../utils/fasy.sh pwasym_stack_cold.txt")  
end program test
