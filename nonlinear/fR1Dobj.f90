module fR1d_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

  COOP_REAL,parameter::coop_fr1d_phi0 = 1.d-6
  COOP_REAL,parameter::coop_fr1d_n = 1.d0
  COOP_REAL,parameter::coop_fr1d_m2cut = 1.d2
  COOP_INT, parameter::coop_fr1d_time_steps = 5
  !!phi = phi0 * (R/R0)**(-n-1)
  !!R = R0 * (phi/phi0)**(-1/(n+1))
  !!R0 = 12 - 9 Omega_m

#define I_NOW this%time(1)
#define I_L this%time(0)
#define I_LL this%time(-1)

  type coop_fr1d_obj
     logical::QS_approx = .true.
     logical::do_GR = .false.
     logical::collapsed = .false.
     COOP_INT::nr, nstep
     !!number of discretized comoving coordinates
     COOP_REAL::rmax, dr, twodr, Omega_m, Omega_Lambda, R0, rho0
     COOP_REAL:: dtau, twodtau
     COOP_REAL,dimension(:),allocatable:: r, dV_out, dV_in, a, Ricci, H, dHdtau, phibg, phibgp
     COOP_REAL,dimension(:,:),allocatable::lnrho, force, u, lnphi, lnphip, lapc
     COOP_INT::time(-1:1)
     COOP_INT, dimension(:,:),allocatable::mask
   contains
     procedure::feedback => coop_fr1d_obj_feedback
     procedure:: free => coop_fr1d_obj_free
     procedure:: init => coop_fr1d_obj_init
     procedure:: rotate_ind => coop_fr1d_obj_rotate_ind
     procedure:: Hofa => coop_fr1d_obj_Hofa
     procedure:: Dofa => coop_fr1d_obj_Dofa
     procedure:: HDofa => coop_fr1d_obj_HDofa
     procedure:: dHdtauofa => coop_fr1d_obj_dHdtauofa
     procedure:: Ricciofa => coop_fr1d_obj_Ricciofa
     procedure:: m2ofphi => coop_fr1d_obj_m2ofphi
     procedure:: deltaRofphi => coop_fr1d_obj_deltaRofphi
     procedure:: phiofdeltaR => coop_fr1d_obj_phiofdeltaR
     procedure:: evolve => coop_fr1d_obj_evolve
     procedure:: evolve_QS => coop_fr1d_obj_evolve_QS
     procedure:: get_force => coop_fr1d_obj_get_force
     procedure:: update_rho => coop_fr1d_obj_update_rho
     procedure:: update_u => coop_fr1d_obj_update_u
     procedure:: update_a => coop_fr1d_obj_update_a
     procedure:: update_phi => coop_fr1d_obj_update_phi
!!     procedure:: get_phi => coop_fr1d_obj_get_phi
  end type coop_fr1d_obj

contains

  subroutine coop_fr1d_obj_rotate_ind(this, backward)
    class(coop_fr1d_obj)::this
    logical,optional::backward
    if(present(backward))then
       if(backward)then
          this%time = mod(this%time + coop_fr1d_time_steps - 1, coop_fr1d_time_steps)
          this%nstep = this%nstep - 1
       else
          this%time = mod(this%time + 1, coop_fr1d_time_steps)
          this%nstep = this%nstep + 1
       endif
    else
       this%time = mod(this%time + 1, coop_fr1d_time_steps)
       this%nstep = this%nstep + 1
    endif
  end subroutine coop_fr1d_obj_rotate_ind

  subroutine coop_fr1d_obj_feedback(this, filename)
    class(coop_fr1d_obj)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file)::fp
    COOP_INT::i
    call fp%open(filename)
    write(fp%unit,"(A2, 4E14.5)") "# ", this%r(0), this%lnrho(0, I_L), this%u(0, I_L), this%lnphi(0, I_L)
    do i=1, this%nr-1
       write(fp%unit,"(7E14.5)") this%r(i), this%lnrho(i, I_L), this%u(i, I_L), this%lnphi(i, I_L), &
            ( (exp(this%lnrho(i, I_L))-this%rho0)/this%a(I_L)**3 - this%deltaRofphi(exp(this%lnphi(i, I_L)), I_L) )*this%a(I_L)**2/3.d0, &
            exp(this%lnphi(i, I_L))*(dot_product(this%lnphi(i-1:i+1, I_L), this%lapc(:, i)) + ((this%lnphi(i+1, I_L)-this%lnphi(i-1, I_L))/this%twodr)**2), &
            exp(this%lnphi(i, I_L))*((this%lnphi(i, I_L)+this%lnphi(i, I_LL)-2.d0*this%lnphi(i, I_L))/this%dtau**2 + (this%lnphi(i, I_NOW)-this%lnphi(i, I_LL))/this%dtau*this%H(I_L) + ((this%lnphi(i, I_NOW)-this%lnphi(i, I_LL))/this%twodtau)**2)
    enddo
    write(fp%unit,"(A2, 4E14.5)") "# ", this%r(this%nr), this%lnrho(this%nr, I_L), this%u(this%nr, I_L), this%lnphi(this%nr, I_L)
    call fp%close()
  end subroutine coop_fr1d_obj_feedback


  subroutine coop_fr1d_obj_free(this)
    class(coop_fr1d_obj)::this
    COOP_DEALLOC(this%r)
    COOP_DEALLOC(this%a)
    COOP_DEALLOC(this%Ricci)
    COOP_DEALLOC(this%H)
    COOP_DEALLOC(this%dHdtau)
    COOP_DEALLOC(this%dV_out)
    COOP_DEALLOC(this%dV_in)
    COOP_DEALLOC(this%lnrho)
    COOP_DEALLOC(this%u)
    COOP_DEALLOC(this%force)
    COOP_DEALLOC(this%mask)
    COOP_DEALLOC(this%lnphi)
    COOP_DEALLOC(this%lnphip)
    COOP_DEALLOC(this%lapc)
  end subroutine coop_fr1d_obj_free

  subroutine coop_fr1d_obj_init(this, Omega_m, nr, rmax, a_ini, delta_ini, r_halo, bw_halo, dtau)
    class(coop_fr1d_obj)::this
    COOP_REAL::omega_m, rmax, rmin, a_ini, delta_ini, r_halo, bw_halo, halfdr, dtau
    COOP_INT:: nr
    COOP_INT::i, j
    call this%free()
    this%rmax = rmax
    this%nr  = nr
    this%nstep = 0
    this%dr = rmax/nr
    this%twodr = this%dr*2.d0
    halfdr = this%dr/2.d0
    this%Omega_m = Omega_m
    this%Omega_Lambda = 1.d0 - Omega_m
    this%R0 = 12.d0 - 9.d0*Omega_m
    this%rho0 = 3.d0*Omega_m

    allocate(this%r(0:nr), this%dV_out(0:nr), this%dV_in(0:nr), this%a(0:coop_fr1d_time_steps-1), this%phibg(0:coop_fr1d_time_steps-1), this%phibgp(0:coop_fr1d_time_steps-1), this%Ricci(0:coop_fr1d_time_steps-1), this%H(0:coop_fr1d_time_steps-1), this%dHdtau(0:coop_fr1d_time_steps-1))
    allocate(this%u(0:nr, 0:coop_fr1d_time_steps-1), this%lnrho(0:nr, 0:coop_fr1d_time_steps-1), this%force(0:nr, 0:coop_fr1d_time_steps-1), this%mask(0:nr, 0:coop_fr1d_time_steps-1))
    allocate(this%lnphi(0:nr, 0:coop_fr1d_time_steps-1), this%lnphip(0:nr, 0:coop_fr1d_time_steps-1), this%lapc(-1:1, 0:nr))

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

    
    this%collapsed = .false.

    this%time = (/ coop_fr1d_time_steps-2, coop_fr1d_time_steps-1, 0 /)
    !!scale factor
    this%a(I_NOW) = a_ini
    this%H(I_NOW) = this%Hofa(a_ini)
    this%dHdtau(I_NOW) = this%dHdtauofa(a_ini)
    this%Ricci(I_NOW) = this%Ricciofa(a_ini)
    this%phibg(I_NOW) = this%phiofdeltaR(0.d0, I_NOW)
    this%phibgp(I_NOW) = 0.d0

    !!log(comoving density)
    this%lnrho(:, I_NOW) = log(this%rho0 *  (1.d0 - delta_ini*(r_halo/rmax)**3 + delta_ini*(tanh((r_halo - this%r)/bw_halo)/2.d0+0.5d0)))
    call this%get_force(I_NOW)
    this%u(:, I_NOW) =  -(1.d0/3.d0*this%HDofa(a_ini)) * (1.d0 - this%rho0/exp(this%lnrho(:, I_NOW)))
    call this%update_phi(I_NOW)
    this%lnphip = 0.d0

    this%dtau = dtau/4.d0
    this%u(:, I_L) = this%u(:, I_NOW)
    this%lnrho(:, I_L) = this%lnrho(:, I_NOW)
    this%force(:, I_L) = this%force(:, I_NOW)
    this%lnphi(:, I_L) = this%lnphi(:, I_NOW)
    this%lnphip(:, I_L) = this%lnphip(:, I_NOW)

    this%a(I_L) = a_ini
    this%H(I_L) = this%H(I_NOW)
    this%dHdtau(I_L) = this%dHdtau(I_NOW)
    this%Ricci(I_L) = this%Ricci(I_NOW)
    this%phibg(I_L) = this%phibg(I_NOW)
    this%phibgp(I_L) = this%phibgp(I_NOW)
    this%twodtau = this%dtau !!first leapfrog => half step
    call this%evolve_QS()

    this%twodtau = dtau/2.d0
    call this%evolve_QS()
    call this%evolve_QS()
    call this%evolve_QS()
    this%u(:, I_L) = this%u(:, 0)
    this%lnrho(:, I_L) = this%lnrho(:, 0)
    this%force(:, I_L) = this%force(:, 0)
    this%lnphi(:, I_L) = this%lnphi(:, 0)
    this%a(I_L) = this%a(0)
    this%H(I_L) = this%H(0)
    this%H(I_L) = this%H(0)
    this%dHdtau(I_L) = this%dHdtau(0)
    this%Ricci(I_L) = this%Ricci(0)
    this%dtau = dtau
    this%twodtau = dtau*2.d0
    call this%evolve_QS()
    do i=1, coop_fr1d_time_steps-3
       call this%evolve_QS()
    enddo
  end subroutine coop_fr1d_obj_init

  subroutine coop_fr1d_obj_evolve_QS(this)
    class(coop_fr1d_obj)::this
    COOP_INT::inow
    call this%rotate_ind()
    call this%update_a(I_NOW)
    call this%update_u(I_NOW)
    call this%update_rho(I_NOW)
    call this%update_phi(I_NOW)
  end subroutine coop_fr1d_obj_evolve_QS


  subroutine coop_fr1d_obj_evolve(this)
    class(coop_fr1d_obj)::this
    call this%evolve_QS()
    if(this%QS_approx .or. this%collapsed)return
  end subroutine coop_fr1d_obj_evolve


  subroutine coop_fr1d_obj_update_phi(this, inow)
    class(coop_fr1d_obj)::this
    COOP_INT,parameter::maxloop = 5
    COOP_INT::inow
    COOP_INT::i, i_l, i_u, nup, loop, il, ill
    COOP_REAL::m2cut, k2max, s, converge, err, fourdr2, lapln, phi, a2by3, m2
    if(this%do_GR) return
    if(this%collapsed)then
       this%lnphi(:, inow) = this%lnphi(:, mod(inow+coop_fr1d_time_steps-1, coop_fr1d_time_steps))
       return
    endif
    il = mod(inow + coop_fr1d_time_steps -1,  coop_fr1d_time_steps)
    ill = mod(il + coop_fr1d_time_steps -1,  coop_fr1d_time_steps)
    k2max = 10.d0/this%dr**2
    m2cut  =  coop_fr1d_m2cut * k2max
    fourdr2 = 4.d0*this%dr**2 
    a2by3 =  this%a(inow)**2/3.d0
    !$omp parallel do private(m2)
    do i=0, this%nr
       this%lnphi(i, inow) = log(this%phiofdeltaR( (exp(this%lnrho(i, inow))-this%rho0)/this%a(inow)**3, inow ))
       m2 = this%m2ofphi(exp(this%lnphi(i, inow)))*a2by3 
       if(m2 .gt. m2cut .or. i.eq. 0  .or. i.eq. this%nr .or. this%nstep .le. 1 )then
          this%mask(i, inow) = 0
       else
          if(m2 .gt. k2max .or. this%QS_approx)then
             this%mask(i, inow) = 1
          else
             this%lnphi(i, inow) = this%lnphi(i, ill) + this%lnphip(i, il)*this%twodtau
             this%lnphip(i, inow) = this%lnphip(i, ill) + ( &
                  (dot_product(this%lnphi(i-1:i+1, il), this%lapc(:, i)) + (this%lnphi(i+1, il)-this%lnphi(i-1, il))**2/fourdr2) &
                  + (this%deltaRofphi(exp(this%lnphi(i, il)), il) - (exp(this%lnrho(i, il))-this%rho0)/this%a(il)**3 ) * this%a(il)**2/3.d0 *exp(-this%lnphi(i, il))  &
                  - 2.d0*this%H(il)*this%lnphip(i, il)-this%lnphip(i, il)**2)*this%twodtau
             this%mask(i, inow) = 2         
          endif
       endif
    enddo
    !$omp end parallel do
    if(this%mask(this%nr-1, inow) .eq. 2)then
       this%mask(this%nr, inow) = 2
       this%lnphi(this%nr, inow) = this%lnphi(this%nr-1, inow)
       this%lnphip(this%nr, inow) = this%lnphip(this%nr-1, inow)
    endif

    if(this%nstep .gt. 1)then
       where(this%mask(:, inow) .eq. 1)
          this%lnphi(:, inow) = this%lnphi(:, il)
       end where
    endif
    i_l = 1
    i_u = this%nr - 1
    nup = count(this%mask(i_l:i_u, inow) .eq. 1)
    if(nup .eq. 0) return
    do while(this%mask(i_l, inow) .ne. 1) 
       i_l = i_l + 1
    enddo
    do while(this%mask(i_u, inow) .ne. 1)
       i_u = i_u - 1
    enddo
    s = 0.d0
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
    this%lnphi(this%nr, inow) =     this%lnphi(this%nr-1, inow)
    if(this%nstep .gt. 1)then
       where(this%mask(:, inow) .ne. 2)
          this%lnphip(:, inow) = (this%lnphi(:, inow) - this%lnphi(:, il))/this%dtau
       end where
    endif
  end subroutine coop_fr1d_obj_update_phi


  subroutine coop_fr1d_obj_update_rho(this, inow)
    class(coop_fr1d_obj)::this
    COOP_INT::il, ill, inow
    if(this%collapsed)then
       this%lnrho(:, inow) = this%lnrho(:, mod(inow+coop_fr1d_time_steps-1, coop_fr1d_time_steps))
       return
    endif
    il = mod(inow + coop_fr1d_time_steps -1,  coop_fr1d_time_steps)
    ill = mod(il + coop_fr1d_time_steps -1,  coop_fr1d_time_steps)
    this%lnrho(1:this%nr-1, inow) = this%lnrho(1:this%nr-1, ill) - (3.d0*this%u(1:this%nr-1, il) + this%r(1:this%nr-1)*(this%u(2:this%nr, il)-this%u(0:this%nr-2, il) + (this%lnrho(2:this%nr, il) - this%lnrho(0:this%nr-2, il))*this%u(1:this%nr-1, il))/(2.d0*this%dr))*this%twodtau
    this%lnrho(0, inow) = this%lnrho(0, ill) - (3.d0*this%twodtau)*this%u(0, il)
    this%lnrho(this%nr, inow) = this%lnrho(this%nr, ill) - (3.d0*this%u(this%nr, il) + this%r(this%nr)*(this%u(this%nr, il)-this%u(this%nr-1, il) + (this%lnrho(this%nr, il) - this%lnrho(this%nr-1, il))*this%u(this%nr, il))/(this%dr))*this%twodtau
    if(any(this%lnrho(:, inow) .gt. log(this%rho0*2.d4)))then
       write(*,"(A, F10.3)") "Halo collapses at z = ", 1.d0/this%a(inow)-1.d0
       this%collapsed = .true.
    endif
    this%lnrho(this%nr-10:this%nr, inow) = sum(this%lnrho(this%nr-10:this%nr-6, inow))/5.d0
    call this%get_force(inow)
  end subroutine coop_fr1d_obj_update_rho


  subroutine coop_fr1d_obj_update_u(this, inow)
    class(coop_fr1d_obj)::this
    COOP_INT::inow, il, ill
    if(this%collapsed)then
       this%u(:, inow) = this%u(:, mod(inow+coop_fr1d_time_steps-1, coop_fr1d_time_steps))
       return
    endif
    il = mod(inow + coop_fr1d_time_steps -1,  coop_fr1d_time_steps)
    ill = mod(il + coop_fr1d_time_steps -1,  coop_fr1d_time_steps)
    this%u(0, inow) = this%u(0, ill) + (this%force(0, il)  - this%H(il)*this%u(0,il) - this%u(0,il)**2)*this%twodtau
    this%u(this%nr,il) = this%u(this%nr,ill) + (this%force(this%nr,il) - this%H(il)*this%u(this%nr,il) - this%u(this%nr,il)**2 &
         - (0.5d0/this%dr)*(this%u(this%nr,il)**2 - this%u(this%nr-1,il)**2)*this%r(this%nr) &
         )*this%twodtau
    if(this%do_GR)then
       this%u(1:this%nr-1, inow) = this%u(1:this%nr-1, ill) + (this%force(1:this%nr-1,il) - this%H(il)*this%u(1:this%nr-1,il) - this%u(1:this%nr-1,il)**2 &
            - (0.25d0/this%dr)*(this%u(2:this%nr,il)**2 - this%u(0:this%nr-2,il)**2)*this%r(1:this%nr-1) &
            )*this%twodtau
    else
       this%u(1:this%nr-1, inow) = this%u(1:this%nr-1, ill) + (this%force(1:this%nr-1,il) - exp(this%lnphi(1:this%nr-1, il))*(this%lnphi(2:this%nr, il)-this%lnphi(0:this%nr-2, il))/(4.d0*this%dr) - this%H(il)*this%u(1:this%nr-1,il) - this%u(1:this%nr-1,il)**2 &
            - (0.25d0/this%dr)*(this%u(2:this%nr,il)**2 - this%u(0:this%nr-2,il)**2)*this%r(1:this%nr-1) &
            )*this%twodtau
    endif
    this%u(this%nr-10:this%nr, inow) = sum(this%u(this%nr-10:this%nr-6, inow))/5.d0
    if(coop_isnan(this%u(:, inow)))then
       write(*,*) this%a(inow)
       stop "u NAN"
    endif

  end subroutine coop_fr1d_obj_update_u


  subroutine coop_fr1d_obj_get_force(this, inow)
    class(coop_fr1d_obj)::this
    COOP_INT::i, inow
    this%force(0, inow) = 0.d0
    do i=1, this%nr
       this%force(i, inow) = this%force(i-1, inow) + (exp(this%lnrho(i-1, inow))-this%rho0)*this%dV_out(i-1) + (exp(this%lnrho(i, inow))-this%rho0) * this%dV_in(i)
    enddo
    this%force(1:this%nr, inow) = - this%force(1:this%nr, inow)/this%r(1:this%nr)**3/this%a(inow)/coop_8pi
    this%force(0, inow) = this%force(1, inow)
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
  
  subroutine coop_fr1d_obj_update_a(this, inow)
    class(coop_fr1d_obj)::this
    COOP_INT::inow, il, ill
    il = mod(inow + coop_fr1d_time_steps-1, coop_fr1d_time_steps)
    ill = mod(il + coop_fr1d_time_steps-1, coop_fr1d_time_steps)
    this%a(inow) = (sqrt(this%a(ill)) + this%H(il)*sqrt(this%a(il))*this%dtau)**2
    this%H(inow) = this%Hofa(this%a(inow))
    this%dHdtau(inow) = this%dHdtauofa(this%a(inow))
    this%Ricci(inow) = 6.d0/this%a(inow)**2*(this%H(inow)**2 + this%dHdtau(inow))
    this%phibg(inow) = this%phiofdeltaR(0.d0, inow)

    if(this%m2ofphi(this%phibg(inow)) .gt. (this%H(inow)/this%a(inow))**2*1.d4)then
       this%phibgp(inow) = (this%phibg(inow) - this%phibg(il))/this%dtau
    else
       this%phibg(inow) = this%phibg(ill) + this%phibgp(il)*this%twodtau
       if(this%phibg(inow) .lt. 0.d0)then
          print*, this%a(inow), this%phibg(ill), this%phibg(il), this%phibg(inow)
          print*, this%a(inow), this%phibgp(ill), this%phibgp(il), this%phibgp(inow)
          stop "background phi < 0 error; please reduce the time step."
       endif
       this%phibgp(inow) = this%phibgp(ill) + (this%deltaRofphi(this%phibg(il), il)*this%a(il)**2/3.d0 - 2.d0*this%H(il)*this%phibgp(il))*this%twodtau
    endif
  end subroutine coop_fr1d_obj_update_a

!!================== model dependent part ==========================

  function coop_fr1d_obj_deltaRofphi(this, phi, inow) result(deltaR)
    class(coop_fr1d_obj)::this   
    COOP_INT::inow
    COOP_REAL::phi, deltaR
    deltaR = this%R0*(phi/coop_fr1d_phi0)**(-1.d0/(coop_fr1d_n+1.d0)) -  this%Ricci(inow)
  end function coop_fr1d_obj_deltaRofphi

!!1/3 d delta R/d phi
  function coop_fr1d_obj_m2ofphi(this, phi) result(m2eff)
    class(coop_fr1d_obj)::this   
    COOP_REAL::phi, m2eff
    m2eff = (1.d0/3.d0/(coop_fr1d_n+1.d0)/coop_fr1d_phi0)*this%R0* (phi/coop_fr1d_phi0)**(-1.d0/(coop_fr1d_n+1.d0)-1.d0)
  end function coop_fr1d_obj_m2ofphi


  function coop_fr1d_obj_phiofdeltaR(this, deltaR, inow) result(phi)
    class(coop_fr1d_obj)::this   
    COOP_INT::inow
    COOP_REAL::phi, deltaR
    phi = coop_fr1d_phi0 * ((deltaR +  this%Ricci(inow))/this%R0)**(-(coop_fr1d_n+1.d0))
  end function coop_fr1d_obj_phiofdeltaR


#undef I_NOW
#undef I_L
#undef I_LL
end module fR1d_mod
