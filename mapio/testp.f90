program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  implicit none
#include "constants.h"
  type(coop_dynamic_array_real)::data_n, data_s
  type(coop_asy)::fig
  COOP_STRING::bin_window_option = "Gaussian"
  COOP_REAL, dimension(:),allocatable::cl_n, cl_s, ells, cl_model
  COOP_INT::i, j, lmax, n, lmin
  COOP_REAL,dimension(:,:),allocatable::weights
  
  call fig%open("maskstrip.txt")
  call fig%init(xlabel = "$\ell$", ylabel = "$\frac{\delta C_b}{C_b}$", caption="$C^{\rm N,S}_b$: hemisphere power spectra", xmin = 10., xmax = 90., ymin = -0.2, ymax = 0.12, doclip = .true.)

  call data_n%load_txt("clsout/clsout_SASYM_0to90.dat")
  n = data_n%nrows
  lmin = 10
  lmax = 90
  allocate(cl_n(n), cl_s(n), ells(n), cl_model(n))
  
  call coop_set_uniform(n, ells, 2.d0, n+1.d0)
  
  allocate(weights(n, n))
  
  weights = 0.d0
  do i = 1, n
     do j=1, n
        weights(j, i) = exp(-((i-j)/10.d0)**2)*(j+0.5d0)/(j+1.d0)/(j+2.d0)*coop_2pi
     enddo
  enddo
  call coop_asy_label(fig, "$C_b(l)\equiv \sum C_{l'} W(l,l')$;\, $W(l, l')\propto (2l'+1) e^{-[(l-l')/10]^2}$", 52., -0.1)     

  do i = 1, n
     weights(:, i) = weights(:, i)/sum(weights(:, i))
  enddo

!!$  call plot(0, 30, "violet", "dotted", 1.)
!!$  call plot(30, 60, "black", "dotted", 1.)
!!$  call plot(60, 90, "blue", "dotted", 1.)
!!$  call plot(0, 90, "red", "solid", 2., legend = "$\frac{C_b^N-C_b^S}{C_b^{\rm \Lambda CDM}$")
!!$
!!$  call plotNS(60, "red", "solid", 2.)
!!$  call plotNS(30, "blue", "dotted", 2.)  
  call plot_dcl("blue", "dashed", 2.)
  call plot_asym("red", "solid", 2.)
  call fig%legend(0.5, 0.89)
  call fig%close()

contains

  subroutine plot_dcl(color, linetype, linewidth)
    COOP_UNKNOWN_STRING::color, linetype
    COOP_SINGLE::linewidth
    call data_n%load_txt("clsout/clsout_NONE.dat")
    do i=1, n
       cl_n(i) = sum(data_n%f(:,3)*weights(:, i))
       cl_model(i) = sum(data_n%f(:, 2)*weights(:, i))       
    enddo
    call fig%curve(ells(lmin-1:lmax-1), cl_n(lmin-1:lmax-1)/cl_model(lmin-1:lmax-1)-1.d0, color=trim(color), linetype=trim(linetype), linewidth = linewidth, legend="$\frac{C_b-C_b^{\rm \Lambda CDM}}{C_b^{\rm \Lambda CDM}}$")
  end subroutine plot_dcl


    subroutine plot_asym(color, linetype, linewidth)
    COOP_UNKNOWN_STRING::color, linetype
    COOP_SINGLE::linewidth
    call data_n%load_txt("clsout/clsout_NASYM_0to90.dat")
    call data_s%load_txt("clsout/clsout_SASYM_0to90.dat")    
    do i=1, n
       cl_n(i) = sum(data_n%f(:,3)*weights(:, i))
       cl_s(i) = sum(data_s%f(:,3)*weights(:, i))       
       cl_model(i) = sum(data_n%f(:, 2)*weights(:, i))       
    enddo
    call fig%curve(ells(lmin-1:lmax-1), (cl_n(lmin-1:lmax-1)-cl_s(lmin-1:lmax-1))/cl_model(lmin-1:lmax-1), color=trim(color), linetype=trim(linetype), linewidth = linewidth, legend="$\frac{C_b^N-C_b^S}{C_b^{\rm \Lambda CDM}}$")
  end subroutine plot_asym



  subroutine plotNS(theta, color, linetype, linewidth)
    COOP_INT::theta
    COOP_UNKNOWN_STRING::color, linetype
    COOP_SINGLE::linewidth
    call data_s%load_txt("clsout/clsout_NASYM_"//COOP_STR_OF(theta)//"to"//COOP_STR_OF(180-theta)//".dat")
    call data_n%load_txt("clsout/clsout_NS_"//COOP_STR_OF(theta)//"to"//COOP_STR_OF(theta)//".dat")
  do i = 1, n
     cl_n(i) = sum(data_n%f(:,3)*weights(:, i))
     cl_s(i) = sum(data_s%f(:,3)*weights(:, i))
     cl_model(i) = sum(data_n%f(:, 2)*weights(:, i))
  enddo
  call fig%curve(ells(lmin-1:lmax-1), (cl_n(lmin-1:lmax-1)- cl_s(lmin-1:lmax-1))/cl_model(lmin-1:lmax-1), color=trim(color), linetype=trim(linetype), linewidth = linewidth, legend="$\theta = "//COOP_STR_OF(theta)//"^\circ$")
    
end subroutine plotNS
  

  subroutine plot(theta1, theta2, color, linetype, linewidth, legend)
    COOP_INT::theta1, theta2
    COOP_UNKNOWN_STRING::color, linetype
    COOP_SINGLE::linewidth
    COOP_UNKNOWN_STRING, optional::legend
    call data_n%load_txt("clsout/clsout_NASYM_"//COOP_STR_OF(theta1)//"to"//COOP_STR_OF(theta2)//".dat")
    call data_s%load_txt("clsout/clsout_SASYM_"//COOP_STR_OF(theta1)//"to"//COOP_STR_OF(theta2)//".dat")
  do i = 1, n
     cl_n(i) = sum(data_n%f(:,3)*weights(:, i))
     cl_s(i) = sum(data_s%f(:,3)*weights(:, i))
     cl_model(i) = sum(data_n%f(:, 2)*weights(:, i))
  enddo
  if(present(legend))then
     call fig%curve(ells(lmin-1:lmax-1), (cl_n(lmin-1:lmax-1)- cl_s(lmin-1:lmax-1))/cl_model(lmin-1:lmax-1), color=trim(color), linetype=trim(linetype), linewidth = linewidth, legend=legend)
  else
     call fig%curve(ells(lmin-1:lmax-1), (cl_n(lmin-1:lmax-1)- cl_s(lmin-1:lmax-1))/cl_model(lmin-1:lmax-1), color=trim(color), linetype=trim(linetype), linewidth = linewidth, legend="$"//COOP_STR_OF(theta1)//"^\circ<\theta<"//COOP_STR_OF(theta2)//"^\circ$")
  end if
    
  end subroutine plot
  
end program test
