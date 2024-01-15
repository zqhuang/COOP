program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  COOP_REAL,parameter::h0 = 0.67d0
  COOP_REAL,parameter::ombh2 = 0.022d0
  COOP_REAL,parameter::omch2 = 0.12d0
  COOP_REAL,parameter::epsilon_s = 0.25d0
  COOP_REAL,parameter::epsilon_inf = 0.2d0
  COOP_REAL,parameter::tau_re = 0.08d0, As = 2.2d-9, ns = 0.96d0
  COOP_REAL:: Qcpl, h, theta0, theta
  type(coop_cosmology_firstorder)::cosmology
  COOP_INT,parameter::lmin = 2, lmax = 2499
  COOP_REAL, dimension(coop_num_Cls, lmin:lmax)::Cls_scalar, Cls_tensor, Cls_lensed
  COOP_REAL, dimension(coop_num_Cls, lmin:lmax)::Cls_scalar_Q0, Cls_tensor_Q0, Cls_lensed_Q0
  COOP_REAL::ells(lmin:lmax)
  COOP_INT l, ik, iQ
  COOP_REAL norm
  type(coop_asy)::fig
  logical::plot_lensed = .true.
  logical::logscale = .true.
  logical::plot_diff = .true.

  norm = 2.72558**2*1.d12
  if(plot_lensed)then
     if(plot_diff)then
        call fig%open("Cl_lensed_Q_diff.txt")
     else
        call fig%open("Cl_lensed_Q_full.txt")        
     endif
  else
     if(plot_diff)then
        call fig%open("Cl_unlensed_Q_diff.txt")
     else
        call fig%open("Cl_unlensed_Q_full.txt")        
     endif
  endif
  if(plot_diff)then
     call fig%init(xlabel = "$\ell$", ylabel = "$C_{\ell}/C_{\ell,Q=0}-1$ ", xlog = logscale  , xmin =1.8, xmax = real(lmax+1), ymin = -0.15, ymax = 0.05, doclip = .true.,  caption="Coupled Dark Energy $Q = \frac{\partial \ln m_{DM}}{\partial \phi}$")      
  else
     call fig%init(xlabel = "$\ell$", ylabel = "$\frac{\ell(\ell+1)C_{\ell}}{2\pi} (\mu K^2)$ ", xlog = logscale , xmin = 1.8, xmax = real(lmax+1),  caption="Coupled Dark Energy $Q = \frac{\partial \ln m_{DM}}{\partial \phi}$")
  endif
  do l=lmin, lmax
     ells(l) = l
  enddo

  do iQ = 0, 3
     Qcpl= iQ*0.1d0
     write(*,*) "doing Q = "//COOP_STR_OF(Qcpl)
     h = h0
     call cosmology%set_standard_cosmology(Omega_b=ombh2/h**2, Omega_c=omch2/h**2, h = h, tau_re = tau_re, As = As, ns = ns, de_Q = Qcpl, de_epsilon_s = epsilon_s, de_epsilon_inf = epsilon_inf)
     theta = cosmology%r_star/cosmology%distlss
     if(iQ .eq. 0) theta0 = theta
     print*, "theta = ", theta
     if(iQ.ne.0)then
        h = h0*(theta0/theta)**4.6
        call cosmology%set_standard_cosmology(Omega_b=ombh2/h**2, Omega_c=omch2/h**2, h = h, tau_re = tau_re, As = As, ns = ns, de_Q = Qcpl, de_epsilon_s = epsilon_s, de_epsilon_inf = epsilon_inf)
        theta = cosmology%r_star/cosmology%distlss        
        print*, "ajusted theta = ", theta        
     endif
     call cosmology%compute_source(0)
     call cosmology%source(0)%get_All_Cls(lmin, lmax, Cls_scalar)
     if(iQ.eq.0)Cls_scalar_Q0 = Cls_scalar  
     if(plot_lensed)then
        call coop_get_lensing_Cls(lmin, lmax, Cls_Scalar, Cls_lensed)
        Cls_lensed = Cls_lensed+Cls_scalar
        if(iQ.eq.0)Cls_lensed_Q0 = Cls_lensed
        if(plot_diff)then
           Cls_lensed = Cls_lensed / Cls_lensed_Q0 - 1.d0
           call fig%curve(ells, Cls_lensed(coop_index_ClTT,:), linewidth=1., color = fig%color(iQ+1), linetype=fig%linetype(iQ+1), legend="$Q = "//COOP_STR_OF(Qcpl)//"$" )           
        else
           call fig%curve(ells, Cls_lensed(coop_index_ClTT,:)*ells*(ells+1.d0)*(norm/coop_2pi), linewidth=1., color = fig%color(iQ+1), linetype=fig%linetype(iQ+1), legend="$Q = "//COOP_STR_OF(Qcpl)//"$" )
        endif
     else
        if(plot_diff)then
           Cls_scalar = Cls_scalar / Cls_scalar_Q0 - 1.d0
           call fig%curve(ells, Cls_scalar(coop_index_ClTT,:), linewidth=1., color = fig%color(iQ+1), linetype=fig%linetype(iQ+1), legend="$Q = "//COOP_STR_OF(Qcpl)//"$")           
        else
           call fig%curve(ells, Cls_scalar(coop_index_ClTT,:)*ells*(ells+1.d0)*(norm/coop_2pi), linewidth=1., color = fig%color(iQ+1), linetype=fig%linetype(iQ+1), legend="$Q = "//COOP_STR_OF(Qcpl)//"$")
        endif
     endif
  enddo

  if(plot_diff)then
     call fig%label("$\epsilon_s = "//COOP_STR_OF(epsilon_s)//"$", 0.1, 0.93, alignment = "right")
     if(plot_lensed)then
        call fig%label("Lensing included", 0.1, 0.86, alignment = "right")
     else
        call fig%label("Lensing NOT included", 0.1, 0.86, alignment = "right")
     endif
     call fig%legend(0.45, 0.4, 1, .false.)     
  else
     call fig%label("$\epsilon_s = "//COOP_STR_OF(epsilon_s)//"$", 0.1, 0.55, alignment = "right")
     if(plot_lensed)then
        call fig%label("Lensing included", 0.1, 0.49, alignment = "right")
     else
        call fig%label("Lensing NOT included", 0.1, 0.49, alignment = "right")
     endif
     call fig%legend(0.1, 0.9, 1, .false.)
  endif
  call fig%close()
  
end program test
