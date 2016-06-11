program collapse
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  !!**********************units*************************************
  !! density: H_0^2 M_p^2 where M_p is the reduced Planck mass
  !! time: 1/H_0
  !! length: c/H_0
  !! field value: M_p
  !!***************************************************************
  COOP_STRING, parameter:: model = "Hu-Sawicki"  !!"negative-powerlaw"
  !----------------------------------------------------------------
  COOP_REAL, parameter::Omega_b = 0.05, Omega_c = 0.25
  COOP_REAL,parameter::Omega_m = Omega_b + Omega_c
  COOP_REAL,parameter::Omega_Lambda = 1. - Omega_m

  !----------------------functions---------------------------------
  type(coop_function)::Vofphi, fofR, intQofphi
  !! Vofphi: V(phi)
  !! fofR:  f(R) 
  !! intQofphi:  \int Q(phi) d\phi where Q(phi) is the coupling (in principle it is phi dependent; for f(R) model Q = 1/sqrt(6) is a constant);

  !------------------parameters used by the models------------------------
  COOP_REAL::n_HuSawicki, c2_HuSawicki, coupling_Q, negative_powerlaw_index
  !----------------- phi and r profiles
  COOP_INT, parameter:: nr = 50
  COOP_REAL:: phi(0:nr+1), r(0:nr+1)
  COOP_REAL:: t
  !!======================Choose model ======================
  select case(trim(model))
  case("Hu-Sawicki")  !!arXiv: 0705.1158
     !!------------ define coupling ---------------------------
     coupling_Q = 1.d0/sqrt(6.d0)
     call intQofphi%init_polynomial( (/ 0.d0, coupling_Q /) )
     !!------------ derive V(phi) from f(R) --------------------------
     Lambda = 3.*Omega_Lambda  !! in unit  H_0^2;
     !!c_1/c_2 = Lambda is fixed; you can change c_2 and n below
     c2_HuSawicki = 12.d0*Omega_Lambda/Omega_m
     n_HuSawicki = 1.  !!Hu and Sawicki used two cases: n=1 and n=4 
     !!define f(R) function; use COOP intrinsic rational function
     call fofR%init_rational( c_up =(/ 2*Lambda /), alpha_up = (/ 0.d0 /), c_down = (/ c2_HuSawicki, 1.d0 /),  alpha_down = (/ n , 0.d0 /) )
     !!convert Jordan-frame f(R) to Einstein-frame V(phi) 
     call coop_convert_fofR_to_Vofphi(fofR, Lambda, Vofphi)  

  case("negative-powerlaw")
     !!------------ define coupling ---------------------------
     coupling_Q = 0.1
     call intQofphi%init_polynomial( (/ 0.d0, coupling_Q /) )
     !!------------ define V(phi) ----------------------------
     negative_powerlaw_index  = -0.1d0
     call Vofphi%init_powerlaw( c = (/ Lambda /), alpha = (/ negative_powerlaw_index /) )
  case default
     write(*,*) "Model: "//trim(model)
     stop "Unknown model."
  end select

  !!================= initialize the cosmology =======================
  call cosmology%init(h=h)
  call cosmology%add_species(coop_baryon(Omega_b))
  call cosmology%add_species(coop_radiation(cosmology%Omega_radiation()))
  call cosmology%add_species(coop_neutrinos_massless(cosmology%Omega_massless_neutrinos_per_species()*(cosmology%Nnu())))  
  call coop_background_add_coupled_DE_with_potential(cosmology, Omega_c = Omega_c, Vofphi = Vofphi, intQofphi = intQofphi,  err = err)
  if(err .ne. 0) stop "Error: bad V(phi) model. cannot initialize the cosmology"

  
  
  

end program collapse
