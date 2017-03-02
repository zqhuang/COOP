module fR1d_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

  COOP_REAL,parameter::coop_fr1d_phi0 = 1.d-6
  COOP_REAL,parameter::coop_fr1d_n = 1.d0
  !!phi = phi0 * (R/R0)**(-n-1)
  !!R = R0 * (phi/phi0)**(-1/(n+1))
  !!R0 = 12 - 9 Omega_m

  type coop_fr1d_obj
     logical::do_GR = .false.
     logical::collapsed = .false.
     COOP_INT::nr, nstep
     !!number of discretized comoving coordinates
     COOP_REAL::rmax, dr,dr2, halfdr, Omega_m, Omega_Lambda, R0, rho0, a, a2, a2by3, a3,  H, dHdtau, Ricci
     COOP_REAL,dimension(:),allocatable:: r, r2, dV_out, dV_in, lnrho, force,  Vslope, u, phi, pi, phieq
     COOP_REAL,dimension(:,:),allocatable::lapc
   contains
     procedure::feedback => coop_fr1d_obj_feedback
     procedure:: free => coop_fr1d_obj_free
     procedure:: init => coop_fr1d_obj_init
     procedure:: Hofa => coop_fr1d_obj_Hofa
     procedure:: Dofa => coop_fr1d_obj_Dofa
     procedure:: HDofa => coop_fr1d_obj_HDofa
     procedure:: dHdtauofa => coop_fr1d_obj_dHdtauofa
     procedure:: Ricciofa => coop_fr1d_obj_RIcciofa
     procedure:: m2ofphi => coop_fr1d_obj_m2ofphi
     procedure:: phiofm2 => coop_fr1d_obj_phiofm2
     procedure:: Rofphi => coop_fr1d_obj_Rofphi
     procedure:: phiofR => coop_fr1d_obj_phiofR
     procedure:: evolve => coop_fr1d_obj_evolve
     procedure:: get_force => coop_fr1d_obj_get_force
     procedure:: update_rho => coop_fr1d_obj_update_rho
     procedure:: update_u => coop_fr1d_obj_update_u
     procedure:: update_a => coop_fr1d_obj_update_a
     procedure:: update_phi => coop_fr1d_obj_update_phi
  end type coop_fr1d_obj

contains


  subroutine coop_fr1d_obj_feedback(this, filename)
    class(coop_fr1d_obj)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file)::fp
    COOP_INT::i
    call fp%open(filename)
    do i= 0, this%nr
       write(fp%unit,"(5E14.5)") this%r(i), this%lnrho(i), this%u(i), this%phi(i), this%pi(i)
    enddo
    call fp%close()
  end subroutine coop_fr1d_obj_feedback


  subroutine coop_fr1d_obj_free(this)
    class(coop_fr1d_obj)::this
    COOP_DEALLOC(this%r)
    COOP_DEALLOC(this%r2)
    COOP_DEALLOC(this%dV_out)
    COOP_DEALLOC(this%dV_in)
    COOP_DEALLOC(this%lapc)
    COOP_DEALLOC(this%lnrho)
    COOP_DEALLOC(this%force)
    COOP_DEALLOC(this%Vslope)
    COOP_DEALLOC(this%u)
    COOP_DEALLOC(this%phi)
    COOP_DEALLOC(this%phieq)
    COOP_DEALLOC(this%pi)
  end subroutine coop_fr1d_obj_free

  subroutine coop_fr1d_obj_init(this, Omega_m, nr, rmax, a_ini, delta_ini, r_halo, bw_halo)
    class(coop_fr1d_obj)::this
    COOP_REAL::omega_m, rmax, rmin, a_ini, delta_ini, r_halo, bw_halo
    COOP_INT:: nr
    COOP_INT::i, j
    call this%free()
    this%rmax = rmax
    this%nr  = nr
    this%nstep = 0
    this%dr = rmax/nr
    this%dr2 = this%dr**2
    this%halfdr = this%dr/2.d0
    this%Omega_m = Omega_m
    this%Omega_Lambda = 1.d0 - Omega_m
    this%R0 = 12.d0 - 9.d0*Omega_m
    this%rho0 = 3.d0*Omega_m
    allocate(this%r2(0:nr), this%r(0:nr), this%dV_out(0:nr), this%dV_in(0:nr), this%lnrho(0:nr), this%force(0:nr), this%Vslope(0:nr), this%u(0:nr), this%lapc(-1:1, 0:nr), this%phi(0:nr), this%pi(0:nr), this%phieq(0:nr))
    this%r = ( (/ ( i, i = 0, nr ) /) ) * this%dr 
    this%r2 = this%r**2
    !!shell volumes
    this%dV_out = (coop_2pi/3.d0*this%dr)*(this%r**2 + (this%r + this%halfdr)*(2.d0*this%r + this%halfdr))
    this%dV_in(0) =  0.d0
    this%dV_in(1:nr) = (coop_2pi/3.d0*this%dr)*(this%r(1:nr)**2 + (this%r(1:nr)- this%halfdr)*(2.d0*this%r(1:nr)-this%halfdr))

    !!Laplacian operator
    this%lapc(:, 0) = 0.d0
    this%lapc(:, nr) = 0.d0
    this%lapc(-1, 1:nr-1) = ((this%r(1:nr-1)-this%halfdr)/this%dr/this%r(1:nr))**2 
    this%lapc(1, 1:nr-1) = ((this%r(1:nr-1)+this%halfdr)/this%dr/this%r(1:nr))**2  
    this%lapc(0, 1:nr-1) = - this%lapc(-1, 1:nr-1) - this%lapc(1, 1:nr-1)
        
    this%collapsed = .false.

    !!scale factor
    this%a = a_ini
    this%a2 = a_ini**2
    this%a2by3 = this%a2/3.d0
    this%a3 = a_ini**3
    this%H = this%Hofa(a_ini)
    this%dHdtau = this%dHdtauofa(a_ini)
    this%Ricci = this%Ricciofa(a_ini)

    !!log(comoving density)
    this%lnrho  = log(this%rho0 *  (1.d0 - delta_ini*(r_halo/rmax)**3 + delta_ini*(tanh((r_halo - this%r)/bw_halo)/2.d0+0.5d0)))
    call this%get_force()
    this%u =  -(1.d0/3.d0*this%HDofa(a_ini)) * (1.d0 - this%rho0/exp(this%lnrho))
    call this%update_phi()
  end subroutine coop_fr1d_obj_init

  subroutine coop_fr1d_obj_evolve(this, dtau)
    class(coop_fr1d_obj)::this
    COOP_REAL::dtau
    this%nstep = this%nstep + 1
    call this%update_u(dtau)
    call this%update_a(dtau/2.d0)
    call this%update_rho(dtau)
    call this%update_phi(dtau)
    call this%update_a(dtau/2.d0)
  end subroutine coop_fr1d_obj_evolve



  subroutine coop_fr1d_obj_update_phi(this, dtau)
    class(coop_fr1d_obj)::this
    COOP_INT,parameter::NSteps = 20
    COOP_REAL,optional::dtau
    COOP_REAL::maxm2, dt, maxphi, minphi, Rmin
    COOP_INT::i , it
    if(this%collapsed)return
    if(this%do_GR)return
    !$omp parallel do
    do i=0, this%nr
       this%Vslope(i) = (exp(this%lnrho(i))-this%rho0)/this%a3 + this%Ricci
       this%phieq(i) = this%phiofR(max(this%Vslope(i), 1.d-99))
    enddo
    !$omp end parallel do
    if(present(dtau) .and. this%nstep .gt. 1)then

       dt = dtau/Nsteps
       maxm2 = 0.1d0/dt**2 - 1.d0/this%dr2
       minphi = this%phiofm2(maxm2/this%a2)
       Rmin =this%Rofphi(minphi)

       this%phi(0)=this%phieq(0)
       this%phi(this%nr)=this%phieq(this%nr)
       this%pi(0) = 0.d0
       this%pi(this%nr) = 0.d0
       !$omp parallel do
       do i=1, this%nr-1
          if(this%phieq(i) .lt. minphi)then
             this%phi(i) = this%phieq(i)
             this%pi(i) = 0.d0
          endif
       enddo
       do it = 1, Nsteps
          !$omp parallel do
          do i=1, this%nr-1
             if(this%phieq(i) .ge. minphi)then
                if(this%phi(i) .gt. minphi)then
                   this%pi(i) = this%pi(i) + (dot_product(this%lapc(:, i), this%phi(i-1:i+1)) - (this%Vslope(i) - this%Rofphi(this%phi(i)))*this%a2by3 )*dt
                else
                   this%pi(i) = this%pi(i) + (dot_product(this%lapc(:, i), this%phi(i-1:i+1)) - ((this%Vslope(i) - Rmin)*this%a2by3  + maxm2*(this%phi(i) - minphi)))*dt
                endif
             endif
          enddo
          !$omp end parallel do
          this%phi(1:this%nr-1) = this%phi(1:this%nr-1) + this%pi(1:this%nr-1)*dt
       enddo

    else
       this%phi = this%phieq
       this%pi = 0.d0
    endif
  end subroutine coop_fr1d_obj_update_phi


  subroutine coop_fr1d_obj_update_rho(this, dtau)
    class(coop_fr1d_obj)::this
    COOP_REAL::dtau
    COOP_INT::i
    if(this%collapsed)return
    this%lnrho(1:this%nr-1) = this%lnrho(1:this%nr-1) - (3.d0*this%u(1:this%nr-1) + this%r(1:this%nr-1)*(this%u(2:this%nr)-this%u(0:this%nr-2) + (this%lnrho(2:this%nr) - this%lnrho(0:this%nr-2))*this%u(1:this%nr-1))/(2.d0*this%dr))*dtau
    this%lnrho(0) = this%lnrho(0) - (3.d0*dtau)*this%u(0)
    this%lnrho(this%nr) = this%lnrho(this%nr-1)  !!drho/dr = 0 boundary condition
    if(.not. all(this%lnrho .lt. log(this%rho0*1.d3)))then
       write(*,"(A, F10.3)") "Halo collapses at z = ", 1.d0/this%a-1.d0
       this%collapsed = .true.
    endif

    call this%get_force()
  end subroutine coop_fr1d_obj_update_rho


  subroutine coop_fr1d_obj_update_u(this, dtau)
    class(coop_fr1d_obj)::this
    COOP_REAL::dtau
    COOP_INT::i
    if(this%collapsed)return
    if(this%do_GR)then
       this%u(1:this%nr-1) = this%u(1:this%nr-1) + (this%force(1:this%nr-1) - this%H*this%u(1:this%nr-1) - this%u(1:this%nr-1)**2 &
            - (0.5d0/this%dr)*(this%u(2:this%nr) - this%u(0:this%nr-2))*this%u(1:this%nr-1)*this%r(1:this%nr-1) &
            )*dtau
    else
       this%u(1:this%nr-1) = this%u(1:this%nr-1) + (this%force(1:this%nr-1) - this%H*this%u(1:this%nr-1) - this%u(1:this%nr-1)**2 &
            - (0.5d0/this%dr)*(this%u(2:this%nr) - this%u(0:this%nr-2))*this%u(1:this%nr-1)*this%r(1:this%nr-1) &
            -  (this%phi(2:this%nr)-this%phi(0:this%nr-2))/this%r(1:this%nr-1)/(4.d0*this%dr) &
            )*dtau
    endif
    this%u(0) = this%u(0) + (this%force(0)  - this%H*this%u(0) - this%u(0)**2)*dtau
    this%u(this%nr) = this%u(this%nr-1)  !!du/dr = 0 boundary condition
  end subroutine coop_fr1d_obj_update_u


  subroutine coop_fr1d_obj_get_force(this)
    class(coop_fr1d_obj)::this
    COOP_INT::i
    this%force(0) = 0.d0
    do i=1, this%nr
       this%force(i) = this%force(i-1) + (exp(this%lnrho(i-1))-this%rho0)*this%dV_out(i-1) + (exp(this%lnrho(i))-this%rho0) * this%dV_in(i)
    enddo
    this%force(1:this%nr) = - this%force(1:this%nr)/this%r(1:this%nr)**3/this%a/coop_8pi
    this%force(0) = this%force(1)
  end subroutine coop_fr1d_obj_get_force

!!================= cosmology background =========================

  function coop_fr1d_obj_Hofa(this, a) result(H)
    class(coop_fr1d_obj)::this
    COOP_REAL::H
    COOP_REAL::a
    H = sqrt(this%Omega_m/a + this%Omega_Lambda*a**2)
  end function coop_fr1d_obj_Hofa

  function coop_fr1d_obj_Dofa(this, a) result(D)
    class(coop_fr1d_obj)::this
    COOP_REAL::D
    COOP_REAL::a
    D = coop_Growth_fitting(this%Omega_m, -1.d0, 1.d0/a-1.d0)
  end function coop_fr1d_obj_Dofa


  function coop_fr1d_obj_HDofa(this, a) result(HD)
    class(coop_fr1d_obj)::this
    COOP_REAL::HD
    COOP_REAL::a
    HD = log(coop_Growth_fitting(this%Omega_m, -1.d0, 1.d0/(a*exp(0.01d0))-1.d0)/coop_Growth_fitting(this%Omega_m, -1.d0, 1.d0/(a*exp(-0.01d0))-1.d0))/0.02d0*this%Hofa(a)
  end function coop_fr1d_obj_HDofa


  function coop_fr1d_obj_dHdtauofa(this, a) result(dHdtau)
    class(coop_fr1d_obj)::this
    COOP_REAL::dHdtau
    COOP_REAL::a
    dHdtau =  -this%Omega_m/a/2.d0  + this%Omega_Lambda*a**2
  end function coop_fr1d_obj_dHdtauofa

  function coop_fr1d_obj_Ricciofa(this, a) result(Ricci)
    class(coop_fr1d_obj)::this   
    COOP_REAL::Ricci
    COOP_REAL::a
    Ricci = 6.d0/a**2*(this%Hofa(a)**2 + this%dHdtauofa(a))
  end function coop_fr1d_obj_Ricciofa
  
  subroutine coop_fr1d_obj_update_a(this, dtau)
    class(coop_fr1d_obj)::this
    COOP_REAL::dtau, atmp, asave
    asave = this%a
    atmp = (sqrt(this%a) + this%Hofa(this%a)*sqrt(this%a)*dtau/4.d0)**2
    this%a = (sqrt(this%a) + this%Hofa(atmp)*sqrt(atmp)*dtau/2.d0)**2
    this%a2 = this%a**2
    this%a2by3 = this%a2/3.d0
    this%a3 = this%a**3
    this%H = this%Hofa(this%a)
    this%dHdtau = this%dHdtauofa(this%a)
    this%Ricci = 6.d0/this%a**2*(this%H**2 + this%dHdtau)
    if(.not. this%do_GR)this%pi = this%pi*(asave/this%a)**2
  end subroutine coop_fr1d_obj_update_a

!!================== model dependent part ==========================


  function coop_fr1d_obj_Rofphi(this, phi) result(R)
    class(coop_fr1d_obj)::this   
    COOP_REAL::phi, R
    R = this%R0*(phi/coop_fr1d_phi0)**(-1.d0/(coop_fr1d_n+1.d0))
  end function coop_fr1d_obj_Rofphi

!!- 1/3 d R/d phi
  function coop_fr1d_obj_m2ofphi(this, phi) result(m2eff)
    class(coop_fr1d_obj)::this   
    COOP_REAL::phi, m2eff
    m2eff = (1.d0/3.d0/(coop_fr1d_n+1.d0)/coop_fr1d_phi0)*this%R0* (phi/coop_fr1d_phi0)**(-1.d0/(coop_fr1d_n+1.d0)-1.d0)
  end function coop_fr1d_obj_m2ofphi

  function coop_fr1d_obj_phiofm2(this, m2) result(phi)
    class(coop_fr1d_obj)::this   
    COOP_REAL::phi, m2
    phi = ((3.d0*(coop_fr1d_n+1.d0)*coop_fr1d_phi0/this%R0)*m2)**(-(coop_fr1d_n+1.d0)/(coop_fr1d_n+2.d0))*coop_fr1d_phi0
  end function coop_fr1d_obj_phiofm2

  function coop_fr1d_obj_phiofR(this, R) result(phi)
    class(coop_fr1d_obj)::this   
    COOP_INT::inow
    COOP_REAL::phi, R
    phi = coop_fr1d_phi0 * (R/this%R0)**(-(coop_fr1d_n+1.d0))
  end function coop_fr1d_obj_phiofR


end module fR1d_mod


