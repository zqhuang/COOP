program Test2
#include "constants.h"    
  use coop_wrapper_utils
  implicit none
  type(coop_asy)::fig
  call fig%open("wsample_npow.txt")
  call fig%init(xlabel="$a$", ylabel="$w$", xmin=0., xmax=1., ymin=-1., ymax=-0.52)

  call fig%plot_file(filename = "womkfigures/npow_02_10.txt", interpolate="LinearLinear", xcol = "1", ycol = "2", color = "gray", linetype="solid", linewidth=1.5, legend="exact solution")
  call fig%plot_file(filename = "womkfigures/npow_02_10.txt", interpolate="LinearLinear", xcol = "1", ycol = "3", color = "red", linetype="dotted", linewidth=2.5, legend="HBKM parametrization")

  call fig%plot_file(filename = "womkfigures/npow_05_m10.txt", interpolate="LinearLinear", xcol = "1", ycol = "2", color = "gray", linetype="solid", linewidth=1.5)
  call fig%plot_file(filename = "womkfigures/npow_05_m10.txt", interpolate="LinearLinear", xcol = "1", ycol = "3", color = "red", linetype="dotted", linewidth=2.5)
  

  call fig%plot_file(filename = "womkfigures/npow_10_10.txt", interpolate="LinearLinear", xcol = "1", ycol = "2", color = "gray", linetype="solid", linewidth=1.5)  
  call fig%plot_file(filename = "womkfigures/npow_10_10.txt", interpolate="LinearLinear", xcol = "1", ycol = "3", color = "red", linetype="dotted", linewidth=2.5)  
  
!!$  call fig%plot_file(filename = "womkfigures/npow_05_m10.txt", interpolate="LinearLinear", xcol = "1", ycol = "2", color = "gray", linetype="solid", linewidth=1.5)  
!!$  call fig%plot_file(filename = "womkfigures/npow_05_m10.txt", interpolate="LinearLinear", xcol = "1", ycol = "3", color = "red", linetype="dotted", linewidth=2.5)
!!$
!!$  
!!$  call fig%plot_file(filename = "womkfigures/npow_05_m20.txt", interpolate="LinearLinear", xcol = "1", ycol = "2", color = "gray", linetype="solid", linewidth=1.5)  
!!$  call fig%plot_file(filename = "womkfigures/npow_05_m20.txt", interpolate="LinearLinear", xcol = "1", ycol = "3", color = "red", linetype="dotted", linewidth=2.5)

  call fig%legend(0.3, 0.92, 1, .true.)

  call fig%label("$V(\phi)\propto \phi^{-\sigma}, \Omega_m=0.31, h=0.68$", 0.3, 0.72)


  call fig%label("$\sigma = 1, \Omega_k=0.1$", 0.35, 0.54)

  call fig%label("$\sigma = 0.5, \Omega_k=-0.1$", 0.35, 0.3)

  call fig%label("$\sigma = 0.2, \Omega_k= 0.1$", 0.35, 0.12)        
  
  call fig%close()
  
end program Test2



