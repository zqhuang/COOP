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
  integer,parameter::num = 50
  type(coop_asy)::figure
  real*8,dimension(num)::phic, b
  real*8,dimension(num, num)::zrec
  real,parameter::ymin = 0.02, ymax = 0.13, xmin = -4.33, xmax = 0.5
  integer::i, j
  call coop_set_uniform(num, b, dble(ymin), dble(ymax))
  call coop_set_uniform(num, phic, dble(xmin), dble(xmax))

  call figure%open("zrec.txt")
  call figure%init(xlabel="$\phi_{\rm c,r}$", ylabel="$\Omega_b h^2$", ymin=ymin, ymax = ymax, xmin = xmin, xmax = xmax)

  do i=1,num
     do j=1, num
        phi_cr = phic(i)
        ombh2 = b(j)
        zrec(i, j) = z_recomb()
     enddo
  enddo

  call figure%density(z=zrec, xmin = dble(xmin), xmax=dble(xmax), ymin=dble(ymin), ymax=dble(ymax), label="$z_{\rm rec}$")
  call figure%close()

end program compute_delta
