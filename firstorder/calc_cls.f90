program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  !!----------------------------------------
  !!wave number k, because COOP uses fixed k arrays, the actual k will be the one that is closest to the following number
  COOP_STRING::output_root = "vis"
  logical::do_cmb_lensing = .true.
#if DO_ZETA_TRANS
  logical::zeta_single_slice = .false.
#endif  
  !!cosmological parameters
  COOP_REAL,parameter::ombh2 = 0.02225d0
  COOP_REAL,parameter::omch2 = 0.1198d0  !!0.12 LCDM  
  COOP_REAL,parameter::hubble = 0.6727d0  !!H0/100
  COOP_REAL,parameter::tau_re = 0.079d0  !!optical depth2
  COOP_REAL,parameter::As = 2.206d-9   !!amplitude
  COOP_REAL, parameter::r = 0.d0  !! tensor/scalar ratio
  COOP_REAL, parameter::ns = 0.9645d0   !!tilt
  COOP_REAL, parameter::Omega_b = ombh2/hubble**2
  COOP_REAL, parameter::Omega_c = omch2/hubble**2
  !!for EFT Dark Energy I have assumed massless neutrinos, if you want to compare with CAMB/CLASS you need to set mnu = 0

  !!DE background EOS
  COOP_REAL, parameter::w0 = -1.d0
  COOP_REAL, parameter::wa = 0.d0    
  
#if DO_EFT_DE  
  !!define the alpha parameters
  COOP_REAL, parameter::alpha_M0 = 0.d0
  COOP_REAL, parameter::alpha_T0 = 0.d0
  COOP_REAL, parameter::alpha_B0 = 0.d0
  COOP_REAL, parameter::alpha_K0 = 0.d0
  COOP_REAL, parameter::alpha_H0 = 0.d0
#elif DO_COUPLED_DE
  COOP_REAL, parameter::Q0 = 0.d0
  COOP_REAL, parameter::Qa = 0.d0  
#endif  
  !!----------------------------------------
  !! declare other variables  
  type(coop_cosmology_firstorder)::cosmology
  type(coop_function)::fwp1, fQ, alphaM, alphaB, alphaK, alphaT, alphaH
  COOP_INT, parameter::lmin = 2, lmax = 2500
  COOP_REAL::Cls(coop_num_Cls, lmin:lmax), tensCls(coop_num_Cls, lmin:lmax), lensedCls(coop_num_Cls, lmin:lmax), ells(lmin:lmax)
  COOP_REAL::norm, lnorm
  COOP_INT::l
  logical success
  type(coop_file)::fp
  !!----------------------------------------
  !!main code
  !!----------------------------------------
  !!DE EOS
  call fwp1%init_polynomial( (/ 1.d0+w0+wa, -wa /) )

#if DO_EFT_DE
  write(*,*) "Dark Energy Model = Effective field theory DE"
  !!initialize alpha functions as  alpha_X(a) = alpha_X0 H_0^2/H(a)^2, where H(a) is LCDM Hubble 
  call generate_function(alpha_M0, alphaM)
  call generate_function(alpha_T0, alphaT)
  call generate_function(alpha_H0, alphaH)
  call generate_function(alpha_B0, alphaB)
  call generate_function(alpha_K0, alphaK)

  !!initialize cosmology
  call cosmology%set_EFT_cosmology(Omega_b=Omega_b, Omega_c=Omega_c, h = hubble, tau_re = tau_re, As = As, ns = ns, wp1 = fwp1, alphaM = alphaM, alphaK = alphaK, alphaB= alphaB, alphaH = alphaH, alphaT = alphaT)
#elif DO_COUPLED_DE
  write(*,*) "Dark Energy Model = Coupled CDM-DE"
  call fQ%init_polynomial( (/ Q0+Qa, -Qa /) )
  call cosmology%set_coupled_DE_cosmology(Omega_b=Omega_b, Omega_c=Omega_c, h = hubble, tau_re = tau_re, As = As, ns = ns, fwp1 = fwp1, fQ = fQ)
#else
  write(*,*) "Dark Energy Model = Lambda (DE w = -1 )"
  call cosmology%set_standard_cosmology(Omega_b=Omega_b, Omega_c=Omega_c, h = hubble, tau_re = tau_re, As = As, ns = ns)
#endif
  !!----------------------------------------
#if DO_ZETA_TRANS
  if(zeta_single_slice)coop_zeta_single_slice_chi = cosmology%distlss
#endif  
  call cosmology%init_source(0)
  call cosmology%compute_source(0, success = success)
  if(.not. success) stop "Solution blows up exponentially; Model is ruled out."
  write(*,*) "linear perturbation done: sigma_8  = ", cosmology%sigma_8  
  call cosmology%source(0)%get_all_cls(lmin, lmax, Cls)
  norm =cosmology%Tcmb()**2*1.d12    
  call fp%open(trim(output_root)//"_scalCls.dat","w")
  write(*,*) "scalar Cls are saved  in "//trim(output_root)//"_scalCls.dat"  
  write(fp%unit, "(A8, 5A16)") "# ell ", "   TT  ",  "   EE  ",  "   TE   ", " PhiPhi ", " TPhi "
  do l = lmin, lmax
     lnorm =  l*(l+1.d0)/coop_2pi*norm
     write(fp%unit, "(I8, 6E16.7)") l, Cls(coop_index_ClTT, l)*lnorm,  Cls(coop_index_ClEE, l)*lnorm,  Cls(coop_index_ClTE, l)*lnorm, Cls(coop_index_ClLenLen, l)*norm*(l*(l+1.d0))**2, Cls(coop_index_ClTLen, l)*norm*(l*(l+1.d0))**1.5
  enddo
  call fp%close()
  if(do_cmb_lensing)then
     call coop_get_lensing_Cls(lmin, lmax, Cls, lensedCls)
     lensedCls = lensedCls + Cls
     call fp%open(trim(output_root)//"_lensedCls.dat","w")
     write(*,*) "lensed Cls are saved  in "//trim(output_root)//"_lensedCls.dat"     
     write(fp%unit, "(A8, 5A16)") "# ell ", "   TT  ",  "   EE  ",  " BB  ", "  TE "
     do l=lmin, lmax
        lnorm = l*(l+1.d0)/coop_2pi*norm
        write(fp%unit, "(I5, 20E16.7)") l, lensedCls(coop_index_ClTT, l)*lnorm, lensedCls(coop_index_ClEE, l)*lnorm,  lensedCls(coop_index_ClBB, l)*lnorm,  lensedCls(coop_index_ClTE, l)*lnorm
     enddo
     call fp%close()
  endif
#if DO_ZETA_TRANS
  call fp%open(trim(output_root)//"_zetaCls.dat","w")
  write(*,*) "zeta Cls are saved  in "//trim(output_root)//"_zetaCls.dat"     
  write(fp%unit, "(A8, 5A16)") "# ell ", "   TT  ",  "  EE  ",  "  ZZ ", " TE ", "  TZ ",  "  EZ  "
  do l=lmin, lmax
     ells(l) = l
     lnorm = l*(l+1.d0)/coop_2pi*norm
     if(zeta_single_slice)then
        write(fp%unit, "(I5, 20E16.7)") l, Cls(coop_index_ClTT, l)*lnorm, Cls(coop_index_ClEE, l)*lnorm,  Cls(coop_index_Clzetazeta, l)*lnorm, Cls(coop_index_ClTE, l)*lnorm, Cls(coop_index_ClTzeta, l)*lnorm, Cls(coop_index_ClEzeta, l)*lnorm, cosmology%Clzetazeta_at_R(l, cosmology%distlss)*lnorm        
     else
        write(fp%unit, "(I5, 20E16.7)") l, Cls(coop_index_ClTT, l)*lnorm, Cls(coop_index_ClEE, l)*lnorm,  Cls(coop_index_Clzetazeta, l)*lnorm, Cls(coop_index_ClTE, l)*lnorm, Cls(coop_index_ClTzeta, l)*lnorm, Cls(coop_index_ClEzeta, l)*lnorm
     endif
  enddo
  call fp%close()
#endif  
  if(r .gt. 0.d0)then  !! do tensor modes
     call cosmology%compute_source(2, success = success)
     call cosmology%source(2)%get_all_cls(lmin, lmax, tensCls)
     write(*,*) "tensor Cls are saved  in "//trim(output_root)//"_tensCls.dat"          
     call fp%open(trim(output_root)//"_tensCls.dat","w")
     write(fp%unit, "(A8, 5A16)") "# ell ", "   TT  ",  "   EE  ",  " BB  ", "  TE "
     do l = lmin, lmax
        lnorm =  l*(l+1.d0)/coop_2pi*norm
        write(fp%unit, "(I8, 6E16.7)") l, tensCls(coop_index_ClTT, l)*lnorm,  tensCls(coop_index_ClEE, l)*lnorm,  tensCls(coop_index_ClBB, l)*lnorm,  tensCls(coop_index_ClTE, l)*lnorm
     enddo
     call fp%close()
     lensedCls = lensedCls + tensCls
     call fp%open(trim(output_root)//"_totCls.dat","w")
     write(*,*) "total Cls are saved  in "//trim(output_root)//"_totCls.dat"     
     write(fp%unit, "(A8, 5A16)") "# ell ", "   TT  ",  "   EE  ",  " BB  ", "  TE "
     do l=lmin, lmax
        lnorm = l*(l+1.d0)/coop_2pi*norm
        write(fp%unit, "(I5, 20E16.7)") l, lensedCls(coop_index_ClTT, l)*lnorm, lensedCls(coop_index_ClEE, l)*lnorm,  lensedCls(coop_index_ClBB, l)*lnorm,  lensedCls(coop_index_ClTE, l)*lnorm
     enddo
     call fp%close()
     
  endif

contains

  subroutine generate_function(alpha0, f)
    COOP_REAL::alpha0
    type(coop_function)::f    
    COOP_REAL,parameter::omega_m0= 0.3
    COOP_REAL,parameter::omega_r0 = 8.d-5
    COOP_INT, parameter::n = 8192
    COOP_REAL::a(n), alpha(n)
    COOP_INT::i
    call coop_set_uniform(n, a, log(coop_min_scale_factor), log(coop_scale_factor_today))
    a = exp(a)
    alpha = alpha0/(1.d0-omega_m0-omega_r0 + omega_m0/a**3 + omega_r0/a**4)
    call f%init(n = n, xmin = a(1), xmax = a(n), f = alpha, xlog = .true., ylog = .false.)
  end subroutine generate_function

  subroutine generate_latevis()
    COOP_INT, parameter::n = 1024
    COOP_REAL::chi(n), vis(n), a, chiend
    COOP_INT::i
#ifdef DO_ZETA_TRANS    
    call coop_set_uniform(n, chi, 0.d0, cosmology%tau0)
    chiend = cosmology%tau0 - cosmology%tauofa(1.d0/(1.d0+cosmology%zre+cosmology%deltaz*5.d0))
    do i=1, n
       if(chi(i) .lt. chiend)then
          a = cosmology%aoftau(cosmology%tau0 - chi(i))
          vis(i) = cosmology%visofa(a)
       else
          vis(i) = 0.d0
       endif
    enddo
    vis = vis/(sum(vis)*(chi(2)-chi(1)))
    call coop_zeta_user_specified_weight%init(n = n, xmin = chi(1), xmax = chi(n), f = vis, method = COOP_INTERPOLATE_LINEAR, name="latevis")
#endif    
  end subroutine generate_latevis


  subroutine generate_earlyvis()
    COOP_INT, parameter::n = 1024
    COOP_REAL::chi(n), vis(n), a, chiend
    COOP_INT::i
#ifdef DO_ZETA_TRANS    
    call coop_set_uniform(n, chi, 0.d0, cosmology%tau0)
    chiend = cosmology%tau0 - cosmology%tauofa(1.d0/(1.d0+cosmology%zre+cosmology%deltaz*5.d0))
    do i=1, n
       if(chi(i) .ge. chiend)then
          a = cosmology%aoftau(cosmology%tau0 - chi(i))
          vis(i) = cosmology%visofa(a)
       else
          vis(i) = 0.d0
       endif
    enddo
    vis = vis/(sum(vis)*(chi(2)-chi(1)))
    call coop_zeta_user_specified_weight%init(n = n, xmin = chi(1), xmax = chi(n), f = vis, method = COOP_INTERPOLATE_LINEAR, name="earlyvis")
#endif    
  end subroutine generate_earlyvis  
  
end program test
