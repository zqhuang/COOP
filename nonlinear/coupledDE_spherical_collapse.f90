program collapse
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  !!**************************************************************
  !!**********Spherical collapse in Einstein Frame ****************
  !!***************************************************************
  !!**********************units************************************
  !! density: H_0^2 M_p^2 where M_p is the reduced Planck mass
  !! time: 1/H_0
  !! length: c/H_0
  !! field value: M_p
  !!***************************************************************
  COOP_STRING, parameter:: model = "Hu-Sawicki"  !!"negative-powerlaw"
  COOP_REAL,parameter::z_init = 3.
  COOP_REAL,parameter::a_init = 1.d0/(1.d0+z_init)
  COOP_REAL,parameter::overdensity_init = 0.01
  COOP_REAL,parameter::radius_init = 0.01    !!unit c/H_0 = 3000/h Mpc
  !----------------------------------------------------------------
  !!for simplicity I am considering a model without baryon and radiation/neutrinos
  COOP_REAL, parameter::Omega_m = 0.3d0    !!Omega_m is defined as lim_{a->0} rho_m a^3 / rho_critical
  COOP_REAL,parameter::Omega_Lambda = 1. - Omega_m
  COOP_REAL, parameter::Lambda = 3.*Omega_Lambda
  type(coop_cosmology_background)::cosmology
  COOP_INT::err
  !----------------------functions---------------------------------
  type(coop_function)::Vofphi, fofR, intQofphi
  !! Vofphi: V(phi)
  !! fofR:  f(R) 
  !! intQofphi:  \int Q(phi) d\phi where Q(phi) is the coupling (in principle it is phi dependent; for f(R) model Q = 1/sqrt(6) is a constant);

  !------------------parameters used by the models------------------------
  COOP_REAL::n_HuSawicki, c2_HuSawicki, coupling_Q, negative_powerlaw_index
  !----------------- phi and r profiles
  COOP_INT, parameter:: nr = 50, nt = 200, n_buffer = 30
  COOP_REAL:: phi(0:nr+n_buffer, 0:nt), phidot(0:nr+n_buffer, 0:nt), r(0:nr+n_buffer, 0:nt), v(0:nr+n_buffer, 0:nt), t(0:nt), z(0:nt), a(0:nt)
  COOP_REAL:: t_end, t_start, dr_init, rhom_init, phi_init, phidot_init, H_init, dt
  !--------------- other variables ---------------------
  COOP_INT::i, j
  !!======================Choose model ======================
  select case(trim(model))
  case("Hu-Sawicki")  !!arXiv: 0705.1158
     !!------------ define coupling ---------------------------
     coupling_Q = 1.d0/sqrt(6.d0)
     call intQofphi%init_polynomial( (/ 0.d0, coupling_Q /) )
     !!------------ derive V(phi) from f(R) --------------------------
     !!c_1/c_2 = Lambda is fixed; you can change c_2 and n below
     c2_HuSawicki = 12.d0*Omega_Lambda/Omega_m
     n_HuSawicki = 1.  !!Hu and Sawicki used two cases: n=1 and n=4 
     !!define f(R) function; use COOP intrinsic rational function
     call fofR%init_rational( c_up =(/ 2*Lambda /), alpha_up = (/ 0.d0 /), c_down = (/ c2_HuSawicki, 1.d0 /),  alpha_down = (/ n_HuSawicki , 0.d0 /) )
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
  call cosmology%init(h=0.7d0)
   call coop_background_add_coupled_DE_with_potential(cosmology, Omega_c = Omega_m, Vofphi = Vofphi, intQofphi = intQofphi,  err = err)
  if(err .ne. 0) stop "Error: bad V(phi) model. cannot initialize the cosmology"
  !!-----------------set time array --------------------------------
  t_start = cosmology%time( a_init )
  t_end = cosmology%time( 1.d0 )
  call coop_set_uniform(nt+1, t, t_start, t_end)
  dt = t(2)-t(1)
  do i = 0, nt
     a(i) = cosmology%aoft(t(i))
     z(i) = 1./a(i) - 1.
  enddo
  phi_init = O0_DE(cosmology)%cplde_phi_lna%eval(log(a_init)) 
  rhom_init = (3.d0*Omega_m)/a_init**3*exp(intQofphi%eval(phi_init))

  !!these are measured at t+dt/2
  H_init = cosmology%Hratio(cosmology%aoft(t(0)+dt/2.d0))
  phidot_init = O0_DE(cosmology)%cplde_phi_prime_lna%eval(log(cosmology%aoft(t(0)+dt/2.d0)))*H_init
  !!================set initial conditions ================================
  !coordinates of test particles; initially uniformly distributed from 0 to radius_init 
  call coop_set_uniform(nr+1, r(0:nr,0), 0.d0, radius_init)
  dr_init = r(1,0) - r(0,0)
  do i= nr+1, nr+n_buffer
     r(i,0) = r(i-1,0)+dr_init
  enddo
  !!create the overdensed region 
  r(0:nr,0) = r(0:nr,0)/(1.d0+overdensity_init)**(1./3.)

  !velocity; leapfrog algorithm: evaluate v at t+dt/2
  v(:,0) = r(:,0) * (cosmology%aoft(t(0)+dt/2.d0)/a(0)) * H_init

  !scalar field
  phi(:,0) = phi_init
  phidot(:, 0) = phidot_init

  print*, phi_init
  do i=0, nt-1
     call evolve(i)
  enddo
  

contains

  !!evolve from step i to i+1
  subroutine evolve(i)
    COOP_INT::i
    r(:, i+1) = r(:, i) + v(:, i)*dt
    phi(:, i+1) = phi(:, i) + phidot(:, i) * dt
!    v(:, i+1) = v(:, i) - 
!    phi(:, i) = phi(:, 
  end subroutine evolve
  

end program collapse
