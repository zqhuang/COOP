program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  !!This program calculates the matter power spectrum
  !!In modified gravity model, matter power spectrum can be ambiguous.
  !!Here I define it as P(k) = k^4 |Phi_k|^2, where Phi is the Newtonian potential in longitudinal gauge.
  !!The output is dimensionless  P(k) in unit of h^{-3} Mpc^3 
  COOP_STRING, parameter::output = "test_matterpower.dat"
  !!----------------------------------------
  !!cosmological parameters
  COOP_REAL,parameter::ombh2 = 0.022d0
  COOP_REAL,parameter::omch2 = 0.12d0  !!0.12 LCDM  
  COOP_REAL,parameter::hubble = 0.676d0  !!H0/100
  COOP_REAL,parameter::tau_re = 0.08d0  !!optical depth2
  COOP_REAL,parameter::As = 2.22d-9   !!amplitude
  COOP_REAL, parameter::ns = 0.96d0   !!tilt
  COOP_REAL, parameter::Omega_b = ombh2/hubble**2
  COOP_REAL, parameter::Omega_c = omch2/hubble**2

  COOP_INT, parameter::nk = 256
  COOP_REAL k(nk), matterPk(nk), khMpc(nk)

  
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
  logical::success
  COOP_INT::ik
  type(coop_file)::fp
  !!----------------------------------------
  !!main code

  !!set a higher feedback level
  coop_feedback_level = 3

  !!initialize w functions as a polynomials of a
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
  !!set k/H0  
  call cosmology%init_source(0)
  call cosmology%compute_source(0, success = success)  
  if(.not. success)write(*,*) "Solution blows up exponentially. Model is ruled out."
  write(*,*) "linear perturbation done: sigma_8  = ", cosmology%sigma_8  
  call coop_set_uniform(nk, k, 0.4d0, 2.d3, logscale = .true.)
  call cosmology%get_Matter_power(z=0.d0, nk = nk, k = k, Pk = matterPk)  !!this returns k^3 |\delta_k|^2 /(2pi^2)
  khMpc = k * cosmology%H0Mpc()/cosmology%h()  !!k/H0 * (H0 * Mpc) / h = k in unit of h Mpc^{-1}
  matterPk = matterPk * (2.d0*coop_pi**2)/khMpc**3  !!obtain |\delta_k|^2 in unit of (Mpc/h)^3
  write(*,*) "Saving the matter power spectrum in "//trim(output)
  call fp%open(trim(output),"w")
  write(fp%unit, "(2A16)") "k * Mpc / h ",  " P(k) * h^3/Mpc^3 "
  do ik=1, nk
     write(fp%unit, "(2E16.7)") khMpc(ik), matterPk(ik)
  enddo
  call fp%close()
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
  
end program test
