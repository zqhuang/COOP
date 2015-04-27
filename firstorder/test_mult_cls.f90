program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  COOP_REAL::Q0 = 0.d0
  COOP_REAL::tracking_n = 0.3d0
  type(coop_cosmology_firstorder)::fod
  COOP_INT,parameter::lmin = 2, lmax = 2000
  COOP_REAL, dimension(coop_num_Cls, lmin:lmax)::Cls_scalar, Cls_tensor, Cls_lensed
  COOP_REAL, dimension(coop_num_Cls, lmin:lmax)::Cls_scalar_Q0, Cls_tensor_Q0, Cls_lensed_Q0
  COOP_REAL::ells(lmin:lmax)
  COOP_INT l, ik
  COOP_REAL norm
  type(coop_file)::fp
  type(coop_asy)::fig
  logical::plot_lensed = .false.
  logical::logscale = .true.
  logical::plot_diff = .true.

  norm = 2.72558**2*1.d12
!$  call fod%set_standard_cosmology(Omega_b=0.0485374d0, Omega_c=0.2585497252d0, h = 0.67766d0, tau_re = 0.08193d0, As = 2.2098d-9, ns = 0.968d0, nrun = 0.d0, r = 0.d0, nt = -0.01d0, YHe = 0.248d0, Nnu = 3.d0, de_Q = Q0, de_tracking_n = tracking_n, de_dlnQdphi = 0.d0, de_dUdphi = 0.d0, de_d2Udphi2 = 0.d0 )
  !!***************************************************
  !! V = V0 / phi^n exp( C  * phi)
  !! n = de_tracking_n;  C = de_dUdphi
  !! V0 is determined by the condition Omega_k =0
  !!***************************************************
  !! Q = Q0 exp( A * phi)
  !! Q0 = de_Q,  A = de_dlnQdphi
  !!***************************************************
!!test energy conservation
!!$  call fod%init_source(0)  
!!$  ik = 1
!!$  do while(fod%source(0)%k(ik).lt. 1.d0)
!!$     ik = ik + 1
!!$  enddo
!!$  call fod%compute_source_k(fod%source(0), ik, do_test_energy_conservation = .true.)
!!$  print*, fod%source(0)%saux(3, ik,  fod%source(0)%ntau),  fod%source(0)%saux(2, ik,  fod%source(0)%ntau)- fod%source(0)%saux(3, ik,  fod%source(0)%ntau)
!!$  
!!$  stop
  
  
  !!compute Cl's
!!$  call fod%compute_source(0)
!!$  call fod%source(0)%get_All_Cls(lmin, lmax, Cls_scalar)
!!$  call fp%open('coupled_DE_scalCls_Q'//trim(coop_num2goodstr(Q0, point="pt"))//'.txt', 'w')
!!$  do l=lmin, lmax
!!$     write(fp%unit, "(I5, 20E16.7)") l, Cls_scalar(:, l)*(l*(l+1.d0)/coop_2pi*norm)
!!$  enddo
!!$  call fp%close()
  if(plot_lensed)then
     call fig%open("Cl_lensed_Q.txt")     
  else
     call fig%open("Cl_Q.txt")
  endif
  if(plot_diff)then
     call fig%init(xlabel = "$\ell$", ylabel = "$\frac{\ell(\ell+1)\delta C_{\ell}}{2\pi} (\mu K^2)$ ", xlog = logscale  , xmin =1.8, xmax = 2200., caption="Coupled Dark Energy $Q = \frac{\partial \ln m_{DM}}{\partial \phi}$")     
  else
     call fig%init(xlabel = "$\ell$", ylabel = "$\frac{\ell(\ell+1)C_{\ell}}{2\pi} (\mu K^2)$ ", xlog = logscale , xmin = 1.8, xmax = 2200., caption="Coupled Dark Energy $Q = \frac{\partial \ln m_{DM}}{\partial \phi}$")
  endif
  do l=lmin, lmax
     ells(l) = l
  enddo
  
!!$  call fod%set_standard_cosmology(Omega_b=0.0485374d0, Omega_c=0.2585497252d0, h = 0.67766d0, tau_re = 0.08193d0, As = 2.2098d-9, ns = 0.968d0, nrun = 0.d0, r = 0.d0, nt = -0.01d0, YHe = 0.248d0, Nnu = 3.d0)
!!$  call fod%compute_source(0)
!!$  call fod%source(0)%get_All_Cls(lmin, lmax, Cls_scalar)
!!$  call fig%curve(ells, Cls_scalar(coop_index_ClTT,:)*ells*(ells+1.d0)*(norm/coop_2pi), linewidth=1., legend="$\Lambda$CDM")

  call fod%set_standard_cosmology(Omega_b=0.0485374d0, Omega_c=0.2585497252d0, h = 0.67766d0, tau_re = 0.08193d0, As = 2.2098d-9, ns = 0.968d0, nrun = 0.d0, r = 0.d0, nt = -0.01d0, YHe = 0.248d0, Nnu = 3.d0, de_Q = 0.d0, de_tracking_n = tracking_n)
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
  
  call fod%set_standard_cosmology(Omega_b=0.0485374d0, Omega_c=0.2585497252d0, h = 0.67766d0, tau_re = 0.08193d0, As = 2.2098d-9, ns = 0.968d0, nrun = 0.d0, r = 0.d0, nt = -0.01d0, YHe = 0.248d0, Nnu = 3.d0, de_Q = 0.1d0, de_tracking_n = tracking_n)
  call fod%compute_source(0)
  call fod%source(0)%get_All_Cls(lmin, lmax, Cls_scalar)
  if(plot_lensed)then
     call coop_get_lensing_Cls(lmin, lmax, Cls_Scalar, Cls_lensed)
     Cls_lensed = Cls_lensed+Cls_scalar
     if(plot_diff)then
        Cls_lensed = Cls_lensed - Cls_lensed_Q0
     endif     
     call fig%curve(ells, Cls_lensed(coop_index_ClTT,:)*ells*(ells+1.d0)*(norm/coop_2pi), linewidth=1., color="blue", linetype="dashed", legend="$Q = 0.1$")
  else
     if(plot_diff)then
        Cls_scalar = Cls_scalar - Cls_scalar_Q0
     endif          
     call fig%curve(ells, Cls_scalar(coop_index_ClTT,:)*ells*(ells+1.d0)*(norm/coop_2pi), linewidth=1., color="blue", linetype="dashed", legend="$Q = 0.1$")
  endif

  call fod%set_standard_cosmology(Omega_b=0.0485374d0, Omega_c=0.2585497252d0, h = 0.67766d0, tau_re = 0.08193d0, As = 2.2098d-9, ns = 0.968d0, nrun = 0.d0, r = 0.d0, nt = -0.01d0, YHe = 0.248d0, Nnu = 3.d0, de_Q = 0.2d0, de_tracking_n = tracking_n)
  call fod%compute_source(0)
  call fod%source(0)%get_All_Cls(lmin, lmax, Cls_scalar)
  if(plot_lensed)then
     call coop_get_lensing_Cls(lmin, lmax, Cls_Scalar, Cls_lensed)
     Cls_lensed = Cls_lensed+Cls_scalar
     if(plot_diff)then
        Cls_lensed = Cls_lensed - Cls_lensed_Q0
     endif          
     call fig%curve(ells, Cls_lensed(coop_index_ClTT,:)*ells*(ells+1.d0)*(norm/coop_2pi), linewidth=1.5,  color="red", linetype="dotted", legend="$Q = 0.2$")
  else
     if(plot_diff)then
        Cls_scalar = Cls_scalar - Cls_scalar_Q0
     endif          
     call fig%curve(ells, Cls_scalar(coop_index_ClTT,:)*ells*(ells+1.d0)*(norm/coop_2pi), linewidth=1.5,  color="red", linetype="dotted", legend="$Q = 0.2$")
  end if

  call fig%label("$V=V_0\phi^{-0.3}$", 0.1, 0.4, alignment = "right")
  if(plot_diff)then
     call fig%legend(0.3, 0.4, 1, .false.)     
  else
     call fig%legend(0.1, 0.9, 1, .false.)
  endif
  call fig%close()
  
end program test
