program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  COOP_REAL::Q0 = 0.d0
  COOP_REAL::alpha = 4.65d0 !!power to adjust h
  COOP_REAL::tracking_n = 0.3d0
  COOP_REAL::h0 = 0.68d0
  COOP_REAL::ombh2 = 0.022d0
  COOP_REAL::omch2 = 0.12d0  
  COOP_REAL::h
  type(coop_cosmology_firstorder)::fod
  COOP_INT,parameter::lmin = 2, lmax = 2499
  COOP_REAL, dimension(coop_num_Cls, lmin:lmax)::Cls_scalar, Cls_tensor, Cls_lensed
  COOP_REAL, dimension(coop_num_Cls, lmin:lmax)::Cls_scalar_Q0, Cls_tensor_Q0, Cls_lensed_Q0
  COOP_REAL::theta, theta0
  COOP_REAL::ells(lmin:lmax)
  COOP_INT l, ik
  COOP_REAL norm
  type(coop_file)::fp
  type(coop_asy)::fig
  logical::plot_lensed = .true.
  logical::logscale = .true.
  logical::plot_diff = .true.

  norm = 2.72558**2*1.d12
  if(plot_lensed)then
     call fig%open("Cl_lensed_Q.txt")     
  else
     call fig%open("Cl_Q.txt")
  endif
  if(plot_diff)then
     call fig%init(xlabel = "$\ell$", ylabel = "$C_{\ell}/C_{\ell,Q=0}-1$ ", xlog = logscale  , xmin =1.8, xmax = real(lmax+1), ymin = -0.05, ymax = 0.05, doclip = .true., caption="Coupled Dark Energy $Q = \frac{\partial \ln m_{DM}}{\partial \phi}$")     
  else
     call fig%init(xlabel = "$\ell$", ylabel = "$\frac{\ell(\ell+1)C_{\ell}}{2\pi} (\mu K^2)$ ", xlog = logscale , xmin = 1.8, xmax = real(lmax+1),  caption="Coupled Dark Energy $Q = \frac{\partial \ln m_{DM}}{\partial \phi}$")
  endif
  do l=lmin, lmax
     ells(l) = l
  enddo
  h = h0
  call fod%set_standard_cosmology(Omega_b=ombh2/h**2, Omega_c=omch2/h**2, h = h, tau_re = 0.08d0, As = 2.2d-9, ns = 0.96d0, de_Q = 0.d0, de_tracking_n = tracking_n)  
  theta0=fod%r_star/fod%distlss
  print*, "theta = ", theta0
  print*, "rho_m at recombination", O0_CDM(fod)%rhoa3_ratio(1.d0/(1.d0+fod%z_star))/O0_CDM(fod)%rhoa3_ratio(1.d-5)  
  call fod%compute_source(0)
  call fod%source(0)%get_All_Cls(lmin, lmax, Cls_scalar)
  Cls_scalar_Q0 = Cls_scalar  
  if(plot_lensed)then
     call coop_get_lensing_Cls(lmin, lmax, Cls_Scalar, Cls_lensed)
     Cls_lensed = Cls_lensed+Cls_scalar
     Cls_lensed_Q0 = Cls_lensed
     if(plot_diff)then
        Cls_lensed = 0.d0
     endif
     call fig%curve(ells, Cls_lensed(coop_index_ClTT,:)*ells*(ells+1.d0)*(norm/coop_2pi), linewidth=1., color = "black", linetype="solid", legend="$Q = 0$")
  else
     if(plot_diff)then
        Cls_scalar = 0.d0
     endif     
     call fig%curve(ells, Cls_scalar(coop_index_ClTT,:)*ells*(ells+1.d0)*(norm/coop_2pi), linewidth=1., color = "black", linetype="solid", legend="$Q = 0$")
  endif

  h = h0
  call fod%set_standard_cosmology(Omega_b=ombh2/h**2, Omega_c=omch2/h**2, h = h, tau_re = 0.08d0, As = 2.2d-9, ns = 0.96d0, de_Q = 0.1d0, de_tracking_n = tracking_n)
  theta=fod%r_star/fod%distlss
  h = (theta0/theta)**alpha*h0
  call fod%set_standard_cosmology(Omega_b=ombh2/h**2, Omega_c=omch2/h**2, h = h, tau_re = 0.08d0, As = 2.2d-9, ns = 0.96d0, de_Q = 0.1d0, de_tracking_n = tracking_n)
  print*, "theta before correction = ", theta
  theta=fod%r_star/fod%distlss    
  print*, "theta after correction= ", theta
  print*, "rho_m at recombination", O0_CDM(fod)%rhoa3_ratio(1.d0/(1.d0+fod%z_star))/O0_CDM(fod)%rhoa3_ratio(1.d-5)
  
  call fod%compute_source(0)
  call fod%source(0)%get_All_Cls(lmin, lmax, Cls_scalar)
  if(plot_lensed)then
     call coop_get_lensing_Cls(lmin, lmax, Cls_Scalar, Cls_lensed)
     Cls_lensed = Cls_lensed+Cls_scalar
     if(plot_diff)then
        Cls_lensed = Cls_lensed/Cls_lensed_Q0-1.d0
        call fig%curve(ells, Cls_lensed(coop_index_ClTT,:), linewidth=1., color="blue", linetype="dashed", legend="$Q = 0.1$")        
     else
        call fig%curve(ells, Cls_lensed(coop_index_ClTT,:)*ells*(ells+1.d0)*(norm/coop_2pi), linewidth=1., color="blue", linetype="dashed", legend="$Q = 0.1$")
     endif
  else
     if(plot_diff)then
        Cls_scalar = Cls_scalar/Cls_scalar_Q0-1.d0
        call fig%curve(ells, Cls_scalar(coop_index_ClTT,:), linewidth=1., color="blue", linetype="dashed", legend="$Q = 0.1$")        
     else         
        call fig%curve(ells, Cls_scalar(coop_index_ClTT,:)*ells*(ells+1.d0)*(norm/coop_2pi), linewidth=1., color="blue", linetype="dashed", legend="$Q = 0.1$")
     endif
  endif

  !!Q = 0.2
  h = h0
  call fod%set_standard_cosmology(Omega_b=ombh2/h**2, Omega_c=omch2/h**2, h = h, tau_re = 0.08d0, As = 2.2d-9, ns = 0.96d0, de_Q = 0.2d0, de_tracking_n = tracking_n)    

  theta= fod%r_star/fod%distlss
  h = (theta0/theta)**alpha*h0
  call fod%set_standard_cosmology(Omega_b=ombh2/h**2, Omega_c=omch2/h**2, h = h, tau_re = 0.08d0, As = 2.2d-9, ns = 0.96d0, de_Q = 0.2d0, de_tracking_n = tracking_n)

  print*, "theta before correction = ", theta  
  theta=fod%r_star/fod%distlss    
  print*, "theta after correction= ", theta
  print*, "rho_m at recombination", O0_CDM(fod)%rhoa3_ratio(1.d0/(1.d0+fod%z_star))/O0_CDM(fod)%rhoa3_ratio(1.d-5)
  
  call fod%compute_source(0)
  call fod%source(0)%get_All_Cls(lmin, lmax, Cls_scalar)
  if(plot_lensed)then
     call coop_get_lensing_Cls(lmin, lmax, Cls_Scalar, Cls_lensed)
     Cls_lensed = Cls_lensed+Cls_scalar
     if(plot_diff)then
        Cls_lensed = Cls_lensed/Cls_lensed_Q0-1.d0
        call fig%curve(ells, Cls_lensed(coop_index_ClTT,:), linewidth=1.5,  color="red", linetype="dotted", legend="$Q = 0.2$")        
     else
        call fig%curve(ells, Cls_lensed(coop_index_ClTT,:)*ells*(ells+1.d0)*(norm/coop_2pi), linewidth=1.5,  color="red", linetype="dotted", legend="$Q = 0.2$")
     endif
  else
     if(plot_diff)then
        Cls_scalar = Cls_scalar/Cls_scalar_Q0-1.d0
        call fig%curve(ells, Cls_scalar(coop_index_ClTT,:), linewidth=1.5,  color="red", linetype="dotted", legend="$Q = 0.2$")        
     else
        call fig%curve(ells, Cls_scalar(coop_index_ClTT,:)*ells*(ells+1.d0)*(norm/coop_2pi), linewidth=1.5,  color="red", linetype="dotted", legend="$Q = 0.2$")
     endif
  end if

  call fig%label("$V=V_0\phi^{-0.3}$", 0.1, 0.4, alignment = "right")
  if(plot_diff)then
     call fig%legend(0.45, 0.4, 1, .false.)     
  else
     call fig%legend(0.1, 0.9, 1, .false.)
  endif
  call fig%close()
  
end program test
