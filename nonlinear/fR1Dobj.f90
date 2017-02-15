module fR1d_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

  COOP_REAL,parameter::coop_fr1d_minq = 1.d-7
  COOP_REAL,parameter::coop_fr1d_phi0 = 1.d-7
  COOP_REAL,parameter::coop_fr1d_n = 1.d0
  COOP_REAL,parameter::coop_fr1d_max_Rbyrho = 3.d0
  !!phi = phi0 * (R/R0)**(-n-1)
  !!R = R0 * (phi/phi0)**(-1/(n+1))
  !!R0 = 12 - 9 Omega_m

  type coop_fr1d_obj
     COOP_INT::nr
     !!number of discretized comoving coordinates
     COOP_INT::ns
     !!number of mass shells
     COOP_REAL::rmax, dr, Omega_m, Omega_Lambda, R0, rho0
     COOP_REAL:: a
     COOP_REAL,dimension(:),allocatable::q, v, m, intm, r, dV_out, dV_in, rho, phi, phip, phiave, gradphi, phic, m2eff
     logical, dimension(:),allocatable::mask
     !!q: comoving positions of mass shells (time dependent)
     !!v = a dq/d tau: rescaled comoving velocities of mass shells (time dependent)
     !!m: masses of shells (fixed)
     !!r: discretized comoving coordinates (fixed)
     !!rho: discretized density
     !!intm: integrated mass
     !!phi: discretized scalar field; phi = 1 - df/dR  = Phi - Psi
     !!phip: a^2 d phi/ d tau
     !!phiave: time averaged phi
   contains
     procedure:: free => coop_fr1d_obj_free
     procedure:: init => coop_fr1d_obj_init
     procedure:: Hofa => coop_fr1d_obj_Hofa
     procedure:: Dofa => coop_fr1d_obj_Dofa
     procedure:: HDofa => coop_fr1d_obj_HDofa
     procedure:: dHdtau => coop_fr1d_obj_dHdtau
     procedure:: Ricci => coop_fr1d_obj_Ricci
     procedure:: m2ofphi => coop_fr1d_obj_m2ofphi
     procedure:: phibg => coop_fr1d_obj_phibg
     procedure:: deltaRofphi => coop_fr1d_obj_deltaRofphi
     procedure:: phiofdeltaR => coop_fr1d_obj_phiofdeltaR
     procedure:: get_rho => coop_fr1d_obj_get_rho
     procedure:: sort_shells => coop_fr1d_obj_sort_shells
     procedure:: integrate_mass => coop_fr1d_obj_integrate_mass
     procedure:: evolve => coop_fr1d_obj_evolve
     procedure:: move_shells => coop_fr1d_obj_move_shells
     procedure:: accelerate_shells => coop_fr1d_obj_accelerate_shells
     procedure:: update_phi => coop_fr1d_obj_update_phi
     procedure:: update_phip => coop_fr1d_obj_update_phip
     procedure:: get_coupled_phi => coop_fr1d_obj_get_coupled_phi
  end type coop_fr1d_obj

contains




  subroutine coop_fr1d_obj_get_rho(this)
    class(coop_fr1d_obj)::this
    COOP_INT::i, j1, j2
    COOP_REAL::fr, frp, lambda, low, up
    call this%sort_shells()
    call this%integrate_mass()
    j1 = 0
    
    do i=1, this%nr
       low = this%r(i) - this%dr/2.d0
       up = this%r(i) + this%dr/2.d0
       do while(this%q(j1+1) .lt. low )
          j1 = j1+ 1
       enddo
       if(j1 .ge. this%ns-10)then
          this%rho(i:) = this%rho(i-1)
          exit
       endif
       j2 = j1 +1 
       do while(this%q(j2).lt. up)
          j2  = j2+1
       enddo
       this%rho(i) = (this%intm(j2)-this%intm(j1))/(this%q(j2)-this%q(j1))/(coop_4pi/3.d0*(this%q(j1)*(this%q(j1)+this%q(j2)) + this%q(j2)**2)) 
    enddo
    this%rho(1) = this%rho(2)
  end subroutine coop_fr1d_obj_get_rho


  subroutine coop_fr1d_obj_sort_shells(this)
    class(coop_fr1d_obj)::this
    COOP_INT::i, j1, j2, j3
    COOP_REAL::vtmp, qtmp
    do i = 2, this%ns
       if( this%q(i) .ge. this%q(i-1) ) cycle
       j1 = 1
       j2 = i-1
       do while(j2-j1 .gt. 1)
          j3 = (j1+j2)/2
          if(this%q(i) .ge. this%q(j3))then
             j1 = j3
          else
             j2 = j3
          endif
       enddo
       vtmp = this%v(i)
       qtmp = this%q(i)
       this%v(j2+1:i) = this%v(j2:i-1)
       this%q(j2+1:i) = this%q(j2:i-1)
       this%v(j2) = vtmp
       this%q(j2) = qtmp
    enddo
  end subroutine coop_fr1d_obj_sort_shells


  subroutine coop_fr1d_obj_integrate_mass(this)
    class(coop_fr1d_obj)::this
    COOP_INT::i
    COOP_REAL::lambda1, lambda2, u, v
    lambda1 = (this%q(1)*2.d0/(this%q(1)+this%q(2)))**3
    this%intm(0) = 0.d0
    this%intm(1) = this%m(1)*lambda1
    do i=2, this%ns-1
       u = (this%q(i+1) - this%q(i-1))/(this%q(i)+this%q(i-1))
       v = (this%q(i) - this%q(i-1))/(this%q(i)+this%q(i-1))
       if(u .gt. 0.d0)then
          lambda2 = (v/u)* ((v**2 + 3.d0*(v+1.d0))/(u**2 + 3.d0*(u+1.d0)))
       else
          lambda2 = 0.5d0
       endif
       this%intm(i) = this%intm(i-1) + this%m(i-1)*(1.d0-lambda1) + this%m(i)*lambda2
       lambda1 = lambda2
    enddo
    lambda2 = 0.5d0 * (1.d0 - (this%ns - 0.5d0)/(2.d0*this%ns*(this%ns-2.d0/3.d0)  + 1.d0))  !!this is the initial weight
    this%intm(this%ns) = this%intm(this%ns-1) + this%m(this%ns-1)*(1.d0-lambda1) + this%m(this%ns)*lambda2
    this%intm(this%ns+1) = this%intm(this%ns) + this%m(this%ns)*(1.d0-lambda2)
  end subroutine coop_fr1d_obj_integrate_mass


  function coop_fr1d_obj_Hofa(this, a) result(H)
    class(coop_fr1d_obj)::this
    COOP_REAL::H
    COOP_REAL, optional::a
    if(present(a))then
       H = sqrt(this%Omega_m/a + this%Omega_Lambda*a**2)
    else
       H = sqrt(this%Omega_m/this%a + this%Omega_Lambda*this%a**2)
    endif
  end function coop_fr1d_obj_Hofa

  function coop_fr1d_obj_Dofa(this, a) result(D)
    class(coop_fr1d_obj)::this
    COOP_REAL::D
    COOP_REAL, optional::a
    if(present(a))then
       D = coop_Growth_fitting(this%Omega_m, -1.d0, 1.d0/a-1.d0)
    else
       D = coop_Growth_fitting(this%Omega_m, -1.d0, 1.d0/this%a-1.d0)
    endif
  end function coop_fr1d_obj_Dofa


  function coop_fr1d_obj_HDofa(this, a) result(HD)
    class(coop_fr1d_obj)::this
    COOP_REAL::HD
    COOP_REAL, optional::a
    if(present(a))then
       HD = log(coop_Growth_fitting(this%Omega_m, -1.d0, 1.d0/(a*exp(0.01d0))-1.d0)/coop_Growth_fitting(this%Omega_m, -1.d0, 1.d0/(a*exp(-0.01d0))-1.d0))/0.02d0*this%Hofa(a)
    else
       HD = log(coop_Growth_fitting(this%Omega_m, -1.d0, 1.d0/(this%a*exp(0.01d0))-1.d0)/coop_Growth_fitting(this%Omega_m, -1.d0, 1.d0/(this%a*exp(-0.01d0))-1.d0))/0.02d0*this%Hofa()
    endif
  end function coop_fr1d_obj_HDofa



  function coop_fr1d_obj_dHdtau(this, a) result(dHdtau)
    class(coop_fr1d_obj)::this
    COOP_REAL::dHdtau
    COOP_REAL, optional::a
    if(present(a))then
       dHdtau =  -this%Omega_m/a/2.d0  + this%Omega_Lambda*a**2
    else
       dHdtau =  -this%Omega_m/this%a/2.d0  + this%Omega_Lambda*this%a**2
    endif
  end function coop_fr1d_obj_dHdtau


  subroutine coop_fr1d_obj_free(this)
    class(coop_fr1d_obj)::this
    COOP_DEALLOC(this%q)
    COOP_DEALLOC(this%m)
    COOP_DEALLOC(this%r)
    COOP_DEALLOC(this%rho)
    COOP_DEALLOC(this%intm)
    COOP_DEALLOC(this%phi)
    COOP_DEALLOC(this%phip)
    COOP_DEALLOC(this%phiave)
    COOP_DEALLOC(this%gradphi)
    COOP_DEALLOC(this%phic)
    COOP_DEALLOC(this%m2eff)
    COOP_DEALLOC(this%mask)
  end subroutine coop_fr1d_obj_free

  subroutine coop_fr1d_obj_init(this, Omega_m, nr, rmax, ns, a_ini, delta_ini, r_halo, bw_halo)
    class(coop_fr1d_obj)::this
    COOP_REAL::omega_m, rmax, a_ini, delta_ini, r_halo, bw_halo
    COOP_INT::ns, nr
    COOP_REAL::dq,  mfac
    COOP_INT::i, j
    call this%free()
    this%rmax = rmax
    this%nr  = nr
    this%ns = ns
    this%dr = rmax/nr
    this%Omega_m = Omega_m
    this%Omega_Lambda = 1.d0 - Omega_m
    this%R0 = 12.d0 - 9.d0*Omega_m
    this%rho0 = 3.d0*Omega_m

    allocate(this%r(nr), this%dV_out(nr), this%dV_in(nr), this%intm(0:ns+1),  this%phi(0:nr+1), this%phiave(0:nr+1), this%phic(nr), this%m2eff(nr), this%mask(nr), this%gradphi(ns), this%phip(nr), this%q(0:ns+1), this%v(ns), this%m(ns), this%rho(nr) )


    this%r(1:this%nr) = ( (/ ( i, i = 1, nr ) /) - 0.5d0 ) * this%dr
    mfac = coop_4pi/3.d0*this%dr/2.d0
    this%dV_in = mfac*(this%r**2 + (this%r-this%dr/2.d0)*(2.d0*this%r-this%dr/2.d0))
    this%dV_out = mfac*(this%r**2 + (this%r+this%dr/2.d0)*(2.d0*this%r+this%dr/2.d0))
    this%a = a_ini
    dq = this%rmax/ns
    this%q(1:ns) = ( (/ ( i, i = 1, ns ) /) - 0.5d0 ) * dq
    this%q(0) = 0.d0
    this%q(ns+1) = rmax


    !set mass
    this%m(1:ns) = (coop_4pi * Omega_m * dq) * (3.d0*this%q(1:ns)**2 + dq**2/4.d0)

    mfac = (1.d0+delta_ini+17.d0/21.d0*delta_ini**2)

    this%v =  ((1.d0-1.d0/mfac)*(a_ini/3.d0)*this%HDofa()) * (tanh((this%q(1:ns)-r_halo)/bw_halo)/2.d0-0.5d0)*this%q(1:ns)
    
    this%m = this%m * (mfac + (1.d0-mfac)*(0.5d0+tanh((this%q(1:ns)-r_halo)/bw_halo)/2.d0))
    call this%get_rho()

    call this%get_coupled_phi()
    this%phiave = this%phi
  end subroutine coop_fr1d_obj_init

  subroutine coop_fr1d_obj_move_shells(this, dlna)
    class(coop_fr1d_obj)::this
    COOP_REAL::dlna
    COOP_INT::i
    this%q(1:this%ns) = max(this%q(1:this%ns) + this%v(1:this%ns)*(dlna/this%Hofa()/this%a), coop_fr1d_minq)
  end subroutine coop_fr1d_obj_move_shells

  subroutine coop_fr1d_obj_accelerate_shells(this, dlna)
    class(coop_fr1d_obj)::this
    COOP_REAL::dlna, dtau, fac
    COOP_INT::i, j
    fac = this%a/this%dr/2.d0
    !$omp parallel do private(j)
    do i=1, this%ns
       j = floor(this%q(i)/this%dr+0.5d0)
       this%gradphi(i) = (this%phi(j+1)-this%phi(j))*fac
    enddo
    !$omp end parallel do
    dtau = dlna/this%Hofa()
    where(this%q(1:this%ns) .gt. coop_fr1d_minq)
       this%v(1:this%ns) = this%v(1:this%ns) - ( this%intm(1:this%ns)/this%q(1:this%ns)**2/coop_8pi &
            - (this%Omega_m/2.d0)*this%q(1:this%ns) &
            + this%gradphi(1:this%ns) )*dtau 
    elsewhere
       this%v(1:this%ns) = 0.d0
    end where
  end subroutine coop_fr1d_obj_accelerate_shells

  subroutine coop_fr1d_obj_update_phi(this, dlna)
    class(coop_fr1d_obj)::this
    COOP_REAL::dlna, dtau, phibg, phic
    COOP_INT::i
    dtau = dlna/this%Hofa()
    phibg = this%phibg()
    !$omp parallel do private(phic)
    do i=1, this%nr
       phic = this%phiofdeltaR( (this%rho(i)-this%rho0)/this%a**3*coop_fr1d_max_Rbyrho )
       this%phi(i) = max(this%phi(i) + this%phip(i) * dtau, phic)
    enddo
    !$omp end parallel do
    this%phi(0) = this%phi(1)
    this%phi(this%nr+1) = this%phi(this%nr)
  end subroutine coop_fr1d_obj_update_phi

  subroutine coop_fr1d_obj_update_phip(this, dlna, friction)
    class(coop_fr1d_obj)::this
    COOP_REAL::dlna, a2by3, a3, twoH, dtau, hub, dr2
    COOP_INT::i
    COOP_REAL,optional::friction
    a2by3 = this%a**2/3.d0
    a3 = this%a**3
    hub = this%Hofa()
    dtau = dlna/hub
    twoH = 2.d0*hub
    if(present(friction))twoH = twoH*(1.d0+friction)
    dr2 = this%dr ** 2
    !$omp parallel do
    do i=2, this%nr-1
       this%phip(i) = this%phip(i) + ( &
            a2by3 * (this%deltaRofphi(this%phi(i)) -  (this%rho(i)-this%rho0) / a3 ) &
            + ( this%r(i+1) * this%phi(i+1) + this%r(i-1) * this%phi(i-1) - 2.d0 * this%r(i) * this%phi(i) )/dr2/this%r(i) &
            - twoH*this%phip(i) &
            ) * dtau
    enddo
    !$omp end parallel do
    this%phip(1) = this%phip(1) + ( &
         a2by3 * ( this%deltaRofphi(this%phi(i)) -  (this%rho(1)-this%rho0) / a3 ) &
         - twoH*this%phip(1) &
         ) * dtau
    this%phip(this%nr) = this%phip(this%nr) + ( &
         a2by3 * ( this%deltaRofphi(this%phi(this%nr)) -  (this%rho(this%nr)-this%rho0) / a3  ) &
         - twoH*this%phip(this%nr) &
         ) * dtau

  end subroutine coop_fr1d_obj_update_phip

  subroutine coop_fr1d_obj_evolve(this, dlna)
    COOP_INT,parameter::Nsteps = 5
    class(coop_fr1d_obj)::this
    COOP_REAL::dlna
    COOP_INT::i
    this%a = this%a*exp(dlna/2.d0)
    call this%move_shells(dlna)
    call this%get_rho()
    call this%accelerate_shells(dlna)
    !$omp parallel do
    do i=1, this%nr
       this%phic(i) = this%phiofdeltaR( (this%rho(i)-this%rho0)/this%a**3 )
       this%m2eff(i) = this%m2ofphi(this%phic(i))
    enddo
    !$omp end parallel do
    this%mask = this%m2eff .lt. (0.1d0/this%dr**2)
    if(any(this%mask))then
       where(this%mask)
          this%phiave = 0.d0
       elsewhere
          this%phi(1:this%nr) = this%phic
          this%phiave = this%phic
       end where
       do i=1, Nsteps
          call this%update_phi(dlna/Nsteps)
          call this%update_phip(dlna/Nsteps)
          where(this%mask)
             this%phiave = this%phiave + this%phi
          end where
       enddo
       where(this%mask)
          this%phiave = this%phiave/Nsteps
       end where
    else
       this%phi(1:this%nr) = this%phic
       this%phi(0) = this%phi(1)
       this%phi(this%nr+1) = this%phi(this%nr)
       this%phiave = this%phic
    endif
    this%a = this%a*exp(dlna/2.d0)
  end subroutine coop_fr1d_obj_evolve


  subroutine coop_fr1d_obj_get_coupled_phi(this)
    class(coop_fr1d_obj)::this
    COOP_INT::i
    COOP_REAL::intm,  a2by6, delta_rho, phib
    !$omp parallel do
    do i=1, this%nr
       this%phi(i) = this%phiofdeltaR( (this%rho(i)-this%rho0)/this%a**3 )
    enddo
    !$omp end parallel do
    this%phi(0) = this%phi(1)
    this%phi(this%nr+1) = this%phi(this%nr)
    this%phip = 0.d0
  end subroutine coop_fr1d_obj_get_coupled_phi


  function coop_fr1d_obj_Ricci(this, a) result(Ricci)
    class(coop_fr1d_obj)::this   
    COOP_REAL::Ricci
    COOP_REAL,optional::a
    if(present(a))then
       Ricci = 6.d0/a**2*(this%Hofa(a)**2 + this%dHdtau(a))
    else
       Ricci = 6.d0/this%a**2*(this%Hofa()**2 + this%dHdtau())
    endif
  end function coop_fr1d_obj_Ricci


!!================== model dependent part ==========================

  function coop_fr1d_obj_deltaRofphi(this, phi) result(deltaR)
    class(coop_fr1d_obj)::this   
    COOP_REAL::phi, deltaR
    deltaR = this%R0*(phi/coop_fr1d_phi0)**(-1.d0/(coop_fr1d_n+1.d0)) -  this%Ricci()
  end function coop_fr1d_obj_deltaRofphi

  function coop_fr1d_obj_m2ofphi(this, phi) result(m2eff)
    class(coop_fr1d_obj)::this   
    COOP_REAL::phi, m2eff
    m2eff = (1.d0/3.d0/(coop_fr1d_n+1.d0)/coop_fr1d_phi0)*this%R0* (phi/coop_fr1d_phi0)**(-1.d0/(coop_fr1d_n+1.d0)-1.d0)
  end function coop_fr1d_obj_m2ofphi


  function coop_fr1d_obj_phiofdeltaR(this, deltaR) result(phi)
    class(coop_fr1d_obj)::this   
    COOP_REAL::phi, deltaR
    phi = coop_fr1d_phi0 * ((deltaR +  this%Ricci())/this%R0)**(-(coop_fr1d_n+1.d0))
  end function coop_fr1d_obj_phiofdeltaR


  function coop_fr1d_obj_phibg(this, a) result(phibg)
    class(coop_fr1d_obj)::this   
    COOP_REAL::phibg
    COOP_REAL,optional::a
    if(present(a))then
       phibg = coop_fr1d_phi0 * (this%Ricci(a)/this%R0)**(-(coop_fr1d_n+1.d0))
    else
       phibg = coop_fr1d_phi0 * (this%Ricci()/this%R0)**(-(coop_fr1d_n+1.d0))
    endif
  end function coop_fr1d_obj_phibg


end module fR1d_mod
