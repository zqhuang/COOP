module Ryskin
  use coop_wrapper_utils  
  implicit none

  real*8::Ombh2 = 0.022 !!this is well measured and required by BBN
  real*8::H0 = 70
  real*8::phi_cr = 0.d0
  real*8::eta = 2.5
  
#define OMEGA_R (4.211e-5/(H0/100.)**2)
#define CURRENT_RATIO   (2.75e-8*Ombh2)
#define OMEGA_M ((1.d0-OMEGA_R*(1.-phi_cr))/4.d0)
  
contains

  function rho_m(z)
    real*8 z, rho_m
    rho_m = (3.*H0**2*OMEGA_M)*(1.+z)**0.75
  end function rho_m

  function rho_r(z)
    real*8 z, rho_r
    rho_r = (3.*H0**2*OMEGA_R)*(1.+z)**(4./(1.-phi_cr))
  end function rho_r

  function rho_de(z)
    real*8 z, rho_de
    rho_de = 3.*rho_m(z) - phi_cr * rho_r(z)
  end function rho_de

  function Hubble(z)
    real*8 z, Hubble
    Hubble = sqrt((rho_m(z)+rho_r(z)+rho_de(z))/3.d0)
  end function Hubble

  function z_recomb()
    real*8 z_recomb    
    integer i
    z_recomb = 1.d3
    do i=1,5
       z_recomb = -13.6*1.6e-19/1.38e-23/2.73/(log(CURRENT_RATIO) + (0.75d0-(3+phi_cr)/(1.-phi_cr))*log(1.+z_recomb))/eta-1.d0
    enddo
  end function z_recomb


  function z_eq()
    real*8 z_eq
    z_eq = exp(log(0.25/OMEGA_R)/(4./(1.-phi_cr)-3./4.))-1.
  end function z_eq


  !!y: delta, dot delta
  subroutine delta_eq(n, a, y, yp)
    integer n
    real*8 a, z, aH, Hub, omm, y(n), yp(n)
    z = 1./a-1.d0
    Hub = Hubble(z)
    aH = a*Hub 
    omm = rho_m(z)/(3.*Hub**2)
    yp(1) = y(2)/aH
    yp(2) = (1.5d0*omm*Hub**2*y(1)-2.*Hub*y(2))/aH
  end subroutine delta_eq

  
end module Ryskin

program compute_delta
  use coop_wrapper_utils
  use Ryskin
  implicit none
  integer,parameter::num = 256
  type(coop_asy)::figure
  type(coop_ode)::ode
  real*8,dimension(num)::phic, delta_rat, exact_delta
  real,parameter::ymin = 1., ymax = 6., xmin = -4.33, xmax = 0.5
  integer::i
  real*8::delta_rec, delta_0, aini, Hini, fac
  aini = 1.d-8
  Hini = Hubble(1./aini-1.)
  fac = 0.001
  call coop_set_uniform(num, phic, xmin+0.01d0, xmax-0.01d0)
  call ode%init(n = 2)
  call figure%open("delta.txt")
  call figure%init(xlabel="$\phi_{\rm c,r}$", ylabel="$\delta_0 / \delta_{\rm rec}$", ymin=ymin, ymax = ymax, xmin = xmin, xmax = xmax)

  Ombh2 = 0.022
  eta = 1.
  do i=1,num
     phi_cr = phic(i)
     call ode%set_initial_conditions(xini = aini, yini = (/ 1.d0, Hini*fac /))
     call ode%evolve(delta_eq, 1.d0/(z_recomb()+1.d0))
     delta_rec = ode%y(1)
     call ode%evolve(delta_eq, 1.d0)
     delta_0 = ode%y(1)
     delta_rat(i) = (1.+min(z_recomb(), z_eq()))**0.2049
     exact_delta(i) = delta_0/delta_rec
  enddo
  call figure%plot(x=phic, y=exact_delta, linewidth=1., color="red", linetype="solid", legend = "$\Omega_bh^2 = 0.022,\, \eta = 1$")
  call figure%plot(x=phic, y=delta_rat, linewidth=1.5, color="red", linetype="dotted")

  Ombh2 = 0.12
  eta = 2.5
  do i=1,num
     phi_cr = phic(i)
     call ode%set_initial_conditions(xini = aini, yini = (/ 1.d0,  Hini*fac /))
     call ode%evolve(delta_eq, 1.d0/(z_recomb()+1.d0))
     delta_rec = ode%y(1)
     call ode%evolve(delta_eq, 1.d0)
     delta_0 = ode%y(1)
     delta_rat(i) = (1.+min(z_recomb(), z_eq()))**0.2049
     exact_delta(i) = delta_0/delta_rec
  enddo
  call figure%plot(x=phic, y=exact_delta, linewidth=1., color="orange", linetype="solid", legend = "$\Omega_bh^2=0.12,\, \eta = 2.5$")  
  call figure%plot(x=phic, y=delta_rat, linewidth=1.5, color="orange", linetype="dotted")
  
  Ombh2 = 0.022
  eta = 2.5
  do i=1,num
     phi_cr = phic(i)
     call ode%set_initial_conditions(xini = aini, yini = (/ 1.d0, Hini*fac /))
     call ode%evolve(delta_eq, 1.d0/(z_recomb()+1.d0))
     delta_rec = ode%y(1)
     call ode%evolve(delta_eq, 1.d0)
     delta_0 = ode%y(1)
     delta_rat(i) = (1.+min(z_recomb(), z_eq()))**0.2049
     exact_delta(i) = delta_0/delta_rec
  enddo
  call figure%plot(x=phic, y=exact_delta, linewidth=1., color="blue", linetype="solid", legend = "$\Omega_bh^2 = 0.022,\, \eta = 2.5$")
  call figure%plot(x=phic, y=delta_rat, linewidth=1.5, color="blue", linetype="dotted")


  Ombh2 = 0.022
  eta = 10.
  do i=1,num
     phi_cr = phic(i)
     call ode%set_initial_conditions(xini = aini, yini = (/ 1.d0, Hini*fac /))
     call ode%evolve(delta_eq, 1.d0/(z_recomb()+1.d0))
     delta_rec = ode%y(1)
     call ode%evolve(delta_eq, 1.d0)
     delta_0 = ode%y(1)
     delta_rat(i) = (1.+min(z_recomb(), z_eq()))**0.2049
     exact_delta(i) = delta_0/delta_rec
  enddo
  call figure%plot(x=phic, y=exact_delta, linewidth=1., color="violet", linetype="solid", legend = "$\Omega_bh^2 = 0.022,\, \eta = 10$")
  call figure%plot(x=phic, y=delta_rat, linewidth=1.5, color="violet", linetype="dotted")
  

  call figure%legend(xratio=0.05, yratio=0.3, cols=1, box=.false.)
  call figure%close()

end program compute_delta
