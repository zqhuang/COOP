program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  implicit none
#include "constants.h"
  type(coop_dynamic_array_real)::data_n, data_s
  type(coop_asy)::fig
  COOP_STRING::bin_window_option = "Gaussian"
  COOP_REAL, dimension(:),allocatable::cl_n, cl_s, ells, cl_model, cl_diff_mean, cl_diff_var
  COOP_INT::i, j, lmax, n, lmin
  COOP_REAL,dimension(:,:),allocatable::weights
  
  call fig%open("clpw.txt")


  call data_n%load_txt("clsout/clsout_NONE.dat")
  n = data_n%nrows
  lmin = 10
  lmax = min(n-9, 100)
  call fig%init(xlabel = "$\ell$", ylabel = "$\delta C_b / C_b^{\Lambda CDM}$", caption="hemisphere power difference", xmin = real(lmin), xmax = real(lmax), ymin = -0.2, ymax = 0.2, doclip = .true.)  
  allocate(cl_n(n), cl_s(n), ells(n), cl_model(n), cl_diff_mean(n), cl_diff_var(n))
  
  call coop_set_uniform(n, ells, 2.d0, n+1.d0)
  
  allocate(weights(n, n))
  
  weights = 0.d0
  do i = 1, n
     do j=1, n
        weights(j, i) = exp(-((i-j)/5.d0)**2)*(j+0.5d0)/(j+1.d0)/(j+2.d0)*coop_2pi
     enddo
  enddo
  call coop_asy_label(fig, "$C_b(l)\equiv \sum C_{l'} W(l,l')$;\, $W(l, l')\propto (2l'+1) e^{-[(l-l')/5]^2}$", 52., -0.18)     

  do i = 1, n
     weights(:, i) = weights(:, i)/sum(weights(:, i))
  enddo

  call plot_cosmic_var(57)
  call plot_asym("ASYM", "red", "solid", 2., "max var. asym. S-N")
  call plot_dcl("blue", "solid", 2.)    
  call plot_asym("DP", "green", "dotted", 1.5, "dipole S-N")
  call plot_asym("GP", "orange", "dotted", 1.5, "Galactic S-N")
  call plot_asym_examples()
  call fig%compact_legend(0.05, 0.92, 2)
  call fig%close()

contains

  subroutine plot_cosmic_var(nsim)
    COOP_INT::isim, nsim
    COOP_INT::i
    cl_diff_mean = 0.d0
    cl_diff_var = 0.d0
    do isim=1, nsim
       call data_n%load_txt("clsout/sim"//COOP_STR_OF(isim)//"_NASYM.dat")
       call data_s%load_txt("clsout/sim"//COOP_STR_OF(isim)//"_SASYM.dat")
       do i=1, n
          cl_n(i) = sum(data_n%f(:,3)*weights(:, i))
          cl_s(i) = sum(data_s%f(:,3)*weights(:, i))       
       enddo
       if(isim .eq. 1)then
          do i=1, n
             cl_model(i) = sum(data_n%f(:, 2)*weights(:, i))       
          enddo
       end if
       cl_diff_mean = cl_diff_mean + (cl_n - cl_s)/cl_model
       cl_diff_var = cl_diff_var + ((cl_n - cl_s)/cl_model)**2
    enddo
    cl_diff_mean = cl_diff_mean/nsim
    cl_diff_var = sqrt(cl_diff_var/nsim - cl_diff_mean**2)
    call fig%band(ells(lmin-1:lmax-1), -2.d0*cl_diff_var(lmin-1:lmax-1), 2.d0*cl_diff_var(lmin-1:lmax-1), colorfill = coop_asy_gray_color(0.7), linecolor = "invisible")    
    call fig%band(ells(lmin-1:lmax-1), -cl_diff_var(lmin-1:lmax-1), cl_diff_var(lmin-1:lmax-1), colorfill = coop_asy_gray_color(0.5), linecolor = "invisible")
  end subroutine plot_cosmic_var

  subroutine plot_dcl(color, linetype, linewidth)
    COOP_UNKNOWN_STRING::color, linetype
    COOP_SINGLE::linewidth
    COOP_INT::i
    call data_n%load_txt("clsout/clsout_NONE.dat")
    do i=1, n
       cl_n(i) = sum(data_n%f(:,3)*weights(:, i))
       cl_model(i) = sum(data_n%f(:, 2)*weights(:, i))       
    enddo
    call fig%curve(ells(lmin-1:lmax-1), cl_model(lmin-1:lmax-1)/cl_model(lmin-1:lmax-1)-1.d0, color="darkgray", linetype = "solid", linewidth = 0.5)    
    call fig%curve(ells(lmin-1:lmax-1), cl_n(lmin-1:lmax-1)/cl_model(lmin-1:lmax-1)-1.d0, color=trim(color), linetype=trim(linetype), linewidth = linewidth, legend="fullsky - $\Lambda$CDM")


  end subroutine plot_dcl


  subroutine plot_asym(axis, color, linetype, linewidth, legend)
    COOP_UNKNOWN_STRING::axis, legend
    COOP_UNKNOWN_STRING::color, linetype
    COOP_SINGLE::linewidth
    COOP_INT::i
    call data_n%load_txt("clsout/clsout_N"//trim(axis)//".dat")
    call data_s%load_txt("clsout/clsout_S"//trim(axis)//".dat")    
    do i=1, n
       cl_n(i) = sum(data_n%f(:,3)*weights(:, i))
       cl_s(i) = sum(data_s%f(:,3)*weights(:, i))       
       cl_model(i) = sum(data_n%f(:, 2)*weights(:, i))       
    enddo
    call fig%curve(ells(lmin-1:lmax-1), (cl_n(lmin-1:lmax-1)-cl_s(lmin-1:lmax-1))/cl_model(lmin-1:lmax-1), color=trim(color), linetype=trim(linetype), linewidth = linewidth, legend=trim(legend))
  end subroutine plot_asym


  subroutine plot_asym_examples()
    COOP_INT::i, isim
    do isim = 1, 10
       call data_n%load_txt("clsout/sim"//COOP_STR_OF(isim)//"_NASYM.dat")
       call data_s%load_txt("clsout/sim"//COOP_STR_OF(isim)//"_SASYM.dat")
       do i=1, n
          cl_n(i) = sum(data_n%f(:,3)*weights(:, i))
          cl_s(i) = sum(data_s%f(:,3)*weights(:, i))       
       enddo
       if(isim .eq. 1)then
          do i=1, n
             cl_model(i) = sum(data_n%f(:, 2)*weights(:, i))       
          enddo
       end if
       if(isim.eq.1)then
          call fig%curve(ells(lmin-1:lmax-1), (cl_n(lmin-1:lmax-1)-cl_s(lmin-1:lmax-1))/cl_model(lmin-1:lmax-1), color="darkgray", linetype = "dotted", linewidth = 1., legend = "sim. examples")          
       else
          call fig%curve(ells(lmin-1:lmax-1), (cl_n(lmin-1:lmax-1)-cl_s(lmin-1:lmax-1))/cl_model(lmin-1:lmax-1), color="darkgray", linetype = "dotted", linewidth = 1.)
       endif
    enddo
  end subroutine plot_asym_examples



  subroutine plotNS(theta, color, linetype, linewidth)
    COOP_INT::theta
    COOP_UNKNOWN_STRING::color, linetype
    COOP_SINGLE::linewidth
    COOP_INT::i
    call data_s%load_txt("clsout/clsout_NASYM_"//COOP_STR_OF(theta)//"to"//COOP_STR_OF(180-theta)//".dat")
    call data_n%load_txt("clsout/clsout_NS_"//COOP_STR_OF(theta)//"to"//COOP_STR_OF(theta)//".dat")
  do i = 1, n
     cl_n(i) = sum(data_n%f(:,3)*weights(:, i))
     cl_s(i) = sum(data_s%f(:,3)*weights(:, i))
     cl_model(i) = sum(data_n%f(:, 2)*weights(:, i))
  enddo
  call fig%curve(ells(lmin-1:lmax-1), (cl_n(lmin-1:lmax-1)- cl_s(lmin-1:lmax-1))/cl_model(lmin-1:lmax-1), color=trim(color), linetype=trim(linetype), linewidth = linewidth, legend="$\theta = "//COOP_STR_OF(theta)//"^\circ$")
    
end subroutine plotNS
  

  subroutine plotcap(theta1, theta2, color, linetype, linewidth, legend)
    COOP_INT::theta1, theta2
    COOP_UNKNOWN_STRING::color, linetype
    COOP_SINGLE::linewidth
    COOP_INT::i
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
    
end subroutine plotcap
  
end program test
