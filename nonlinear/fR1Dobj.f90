module fR1d_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

  COOP_REAL,parameter::coop_fr1d_phi0 = 1.d-7
  COOP_REAL,parameter::coop_fr1d_n = 1.d0
  COOP_REAL,parameter::coop_fr1d_m2cut = 1.d2
  !!phi = phi0 * (R/R0)**(-n-1)
  !!R = R0 * (phi/phi0)**(-1/(n+1))
  !!R0 = 12 - 9 Omega_m

#define PHI_NOW this%ind_phi(1)
#define PHI_L this%ind_phi(0)
#define PHI_LL this%ind_phi(-1)
#define U_NOW this%ind_u(1)
#define U_L this%ind_u(0)


  type coop_fr1d_obj
     logical::QS_approx = .true.
     logical::ignore_cameleon_force = .false.
     COOP_INT::nr
     !!number of discretized comoving coordinates
     COOP_REAL::rmax, dr, Omega_m, Omega_Lambda, R0, rho0
     COOP_REAL:: a, alast
     COOP_REAL:: dtau
     COOP_REAL,dimension(:),allocatable:: r, dV_out, dV_in
     COOP_REAL,dimension(:,:),allocatable::lnrho, force, u, lnphi, lapc
     COOP_INT::ind_phi(-1:1), ind_u(0:1)
     logical, dimension(:,:),allocatable::mask
   contains
     procedure:: free => coop_fr1d_obj_free
     procedure:: init => coop_fr1d_obj_init
     procedure:: rotate_ind => coop_fr1d_obj_rotate_ind
     procedure:: Hofa => coop_fr1d_obj_Hofa
     procedure:: Dofa => coop_fr1d_obj_Dofa
     procedure:: HDofa => coop_fr1d_obj_HDofa
     procedure:: dHdtau => coop_fr1d_obj_dHdtau
     procedure:: Ricci => coop_fr1d_obj_Ricci
     procedure:: m2ofphi => coop_fr1d_obj_m2ofphi
     procedure:: phibg => coop_fr1d_obj_phibg
     procedure:: deltaRofphi => coop_fr1d_obj_deltaRofphi
     procedure:: phiofdeltaR => coop_fr1d_obj_phiofdeltaR
     procedure:: evolve => coop_fr1d_obj_evolve
     procedure:: get_force => coop_fr1d_obj_get_force
     procedure:: update_rho => coop_fr1d_obj_update_rho
     procedure:: update_u => coop_fr1d_obj_update_u
     procedure:: update_a => coop_fr1d_obj_update_a
     procedure:: get_QS_phi => coop_fr1d_obj_get_QS_phi
  end type coop_fr1d_obj

contains

  subroutine coop_fr1d_obj_rotate_ind(this)
    class(coop_fr1d_obj)::this
    this%ind_phi = mod(this%ind_phi + 1, 3)
    this%ind_u = mod(this%ind_u + 1, 2)
  end subroutine coop_fr1d_obj_rotate_ind

  subroutine coop_fr1d_obj_free(this)
    class(coop_fr1d_obj)::this
    COOP_DEALLOC(this%r)
    COOP_DEALLOC(this%dV_out)
    COOP_DEALLOC(this%dV_in)
    COOP_DEALLOC(this%lnrho)
    COOP_DEALLOC(this%u)
    COOP_DEALLOC(this%force)
    COOP_DEALLOC(this%mask)
    COOP_DEALLOC(this%lnphi)
    COOP_DEALLOC(this%lapc)
  end subroutine coop_fr1d_obj_free

  subroutine coop_fr1d_obj_init(this, Omega_m, nr, rmax, a_ini, delta_ini, r_halo, bw_halo, dtau)
    class(coop_fr1d_obj)::this
    COOP_REAL::omega_m, rmax, rmin, a_ini, delta_ini, r_halo, bw_halo, halfdr, dtau
    COOP_INT:: nr
    COOP_INT::i, j
    call this%free()
    this%ind_phi = (/ 0, 1, 2 /)
    this%ind_u = (/ 0, 1 /)
    this%rmax = rmax
    this%nr  = nr
    this%dr = rmax/nr
    halfdr = this%dr/2.d0
    this%Omega_m = Omega_m
    this%Omega_Lambda = 1.d0 - Omega_m
    this%R0 = 12.d0 - 9.d0*Omega_m
    this%rho0 = 3.d0*Omega_m

    allocate(this%r(0:nr), this%dV_out(0:nr), this%dV_in(0:nr))
    allocate(this%u(0:nr, 0:1), this%lnrho(0:nr, 0:1), this%force(0:nr, 0:1), this%mask(0:nr, 0:1))
    allocate(this%lnphi(0:nr, 0:2), this%lapc(-1:1, 0:nr))

    this%r = ( (/ ( i, i = 0, nr ) /) ) * this%dr 
    !!shell volumes
    this%dV_out = (coop_2pi/3.d0*this%dr)*(this%r**2 + (this%r + halfdr)*(2.d0*this%r + halfdr))
    this%dV_in(0) =  0.d0
    this%dV_in(1:nr) = (coop_2pi/3.d0*this%dr)*(this%r(1:nr)**2 + (this%r(1:nr)- halfdr)*(2.d0*this%r(1:nr)-halfdr))
    !!Laplacian operator
    this%lapc(:, 0) = 0.d0
    this%lapc(:, nr) = 0.d0
    this%lapc(-1, 1:nr-1) = this%r(0:nr-2)/this%dr**2/this%r(1:nr-1)
    this%lapc(1, 1:nr-1) = this%r(2:nr)/this%dr**2/this%r(1:nr-1)
    this%lapc(0, 1:nr-1) = -2.d0/this%dr**2
    !!scale factor
    this%a = a_ini
    !!time step
    this%dtau = dtau
    !!log(comoving density)
    this%lnrho(:, U_NOW) = log(this%rho0 *  (1.d0+delta_ini*(tanh((r_halo - this%r)/bw_halo)/2.d0+0.5d0)))
    call this%get_force(rho_index = U_NOW, force_index = U_NOW)
    this%u(:, U_NOW) =  -(1.d0/3.d0*this%HDofa()) * (1.d0 - this%rho0/exp(this%lnrho(:, U_NOW)))
    call this%get_QS_phi(rho_index = U_NOW, phi_index = PHI_NOW, mask_index = U_NOW)
    call this%evolve()


  end subroutine coop_fr1d_obj_init

  subroutine coop_fr1d_obj_evolve(this)
    class(coop_fr1d_obj)::this
    call this%rotate_ind()
    call this%update_u(index_from = U_L, index_to = U_NOW)
    call this%update_rho(index_from = U_L, index_to = U_NOW)
    call this%update_a()
    call this%get_QS_phi(rho_index = U_NOW, phi_index = PHI_NOW, mask_index = U_NOW)
  end subroutine coop_fr1d_obj_evolve




  subroutine coop_fr1d_obj_get_QS_phi(this, rho_index, phi_index, mask_index)
    class(coop_fr1d_obj)::this
    COOP_INT,parameter::maxloop = 20
    COOP_INT::phi_index, mask_index, rho_index
    COOP_INT::i, i_l, i_u, nup, loop
    COOP_REAL::m2cut, s, converge, err, fourdr2, lapln, phi, a2by3
    m2cut =  coop_fr1d_m2cut/this%dr**2
    !$omp parallel do
    do i=0, this%nr
       this%lnphi(i, phi_index) = log(this%phiofdeltaR( (exp(this%lnrho(i, rho_index))-this%rho0)/this%a**3 ))
       this%mask(i, mask_index) = this%m2ofphi(exp(this%lnphi(i, phi_index))) .lt. m2cut
    enddo
    !$omp end parallel do
    i_l = 1
    i_u = this%nr - 1
    nup = count(this%mask(i_l:i_u, mask_index))
    if(nup .eq. 0) return
    s = 0.d0
    fourdr2 = 4.d0*this%dr**2 
    a2by3 =  this%a**2/3.d0
    do i = i_l, i_u
#include "updatephi_QS.h"
    enddo
    do i = i_u, i_l, -1
#include "updatephi_QS.h"       
    enddo
    converge = s/nup/2.d0 * 1.d-4
    loop = 0
    do while(s/nup/2.d0 .gt. converge .and. loop .lt. maxloop)
       s = 0.d0
       do i = i_l, i_u
#include "updatephi_QS.h"
       enddo
       do i = i_u, i_l, -1
#include "updatephi_QS.h"       
       enddo
       loop = loop + 1
    enddo
  end subroutine coop_fr1d_obj_get_QS_phi


  subroutine coop_fr1d_obj_update_rho(this, index_from, index_to)
    class(coop_fr1d_obj)::this
    COOP_INT::index_from, index_to
    this%lnrho(1:this%nr-1, index_to) = this%lnrho(1:this%nr-1, index_from) - (3.d0*this%u(1:this%nr-1, index_from) + this%r(1:this%nr-1)*(this%u(2:this%nr, index_from)-this%u(0:this%nr-2, index_from) + (this%lnrho(2:this%nr, index_from) - this%lnrho(0:this%nr-2, index_from))*this%u(1:this%nr-1, index_from))/(2.d0*this%dr))*this%dtau
    this%lnrho(0, index_to) = this%lnrho(0, index_from) - (3.d0*this%dtau)*this%u(0, index_from)
    this%lnrho(this%nr, index_to) = this%lnrho(this%nr, index_from) - (3.d0*this%dtau)*this%u(this%nr, index_from)
    if(any(this%lnrho(:, index_to) .gt. log(this%rho0)+10.d0))then
       write(*,"(A, F10.3)") "Halo collapses at z = ", 1.d0/this%a-1.d0
       stop 
    endif
    call this%get_force(rho_index = index_to, force_index = index_to)
  end subroutine coop_fr1d_obj_update_rho


  subroutine coop_fr1d_obj_update_u(this, index_from, index_to)
    class(coop_fr1d_obj)::this
    COOP_INT::index_from, index_to
    COOP_REAL::H
    H = this%Hofa()
    this%u(0, index_to) = this%u(0, index_from) + (this%force(0, index_from)  - H*this%u(0,index_from) - this%u(0,index_from)**2)*this%dtau
    this%u(this%nr,index_from) = this%u(this%nr,index_from) + (this%force(this%nr,index_from) - H*this%u(this%nr,index_from) - this%u(this%nr,index_from)**2)*this%dtau
    this%u(1:this%nr-1, index_to) = this%u(1:this%nr-1, index_from) + (this%force(1:this%nr-1,index_from) - exp(this%lnphi(1:this%nr-1, index_from))*(this%lnphi(2:this%nr, index_from)-this%lnphi(0:this%nr-2, index_from))/(2.d0*this%dr) - H*this%u(1:this%nr-1,index_from) - this%u(1:this%nr-1,index_from)**2 &
         - (0.25d0/this%dr)*(this%u(2:this%nr,index_from)**2 - this%u(0:this%nr-2,index_from)**2)*this%r(1:this%nr-1) &
         )*this%dtau
    if(coop_isnan(this%u(:, index_to)))then
       write(*,*) this%a
       stop "u NAN"
    endif

  end subroutine coop_fr1d_obj_update_u


  subroutine coop_fr1d_obj_get_force(this, rho_index, force_index)
    class(coop_fr1d_obj)::this
    COOP_INT::rho_index, force_index, i
    this%force(0, force_index) = 0.d0
    do i=1, this%nr
       this%force(i, force_index) = this%force(i-1, force_index) + (exp(this%lnrho(i-1, rho_index))-this%rho0)*this%dV_out(i-1) + (exp(this%lnrho(i, rho_index))-this%rho0) * this%dV_in(i)
    enddo
    this%force(1:this%nr, force_index) = - this%force(1:this%nr, force_index)/this%r(1:this%nr)**3/this%a/coop_8pi
    this%force(0, force_index) = this%force(1, force_index)
  end subroutine coop_fr1d_obj_get_force

!!================= cosmology background =========================

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

  subroutine coop_fr1d_obj_update_a(this)
    class(coop_fr1d_obj)::this
    COOP_REAL::da
    this%alast= this%a
    da = this%Hofa()*this%a*this%dtau
    da = this%Hofa(this%a+da/2.d0)*(this%a+da/2.d0)*this%dtau
    this%a = this%a + this%Hofa(this%a+da/2.d0)*(this%a+da/2.d0)*this%dtau 
  end subroutine coop_fr1d_obj_update_a

!!================== model dependent part ==========================

  function coop_fr1d_obj_deltaRofphi(this, phi) result(deltaR)
    class(coop_fr1d_obj)::this   
    COOP_REAL::phi, deltaR
    deltaR = this%R0*(phi/coop_fr1d_phi0)**(-1.d0/(coop_fr1d_n+1.d0)) -  this%Ricci()
  end function coop_fr1d_obj_deltaRofphi

!!1/3 d delta R/d phi
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

#undef PHI_NOW
#undef PHI_L
#undef PHI_LL
#undef U_NOW
#undef U_L
end module fR1d_mod
