module coop_firstorder_mod
  use coop_wrapper_background
  use coop_recfast_mod
  use coop_pertobj_mod
  implicit none
#include "constants.h"

!private

  public::coop_cosmology_firstorder, coop_cosmology_firstorder_source
  
  COOP_REAL, parameter :: coop_power_lnk_min = log(0.1d0) 
  COOP_REAL, parameter :: coop_power_lnk_max = log(5.d3) 
  COOP_REAL, parameter :: coop_visibility_amin = 1.8d-4
  COOP_REAL, parameter :: coop_initial_condition_epsilon = 1.d-6
  COOP_REAL, parameter :: coop_cosmology_firstorder_ode_accuracy = 1.d-7
  COOP_REAL, parameter :: coop_cosmology_firstorder_tc_cutoff = 0.06d0

  COOP_REAL, dimension(0:2), parameter::coop_source_tau_weight = (/ 0.15d0, 0.15d0, 0.1d0 /)
  COOP_INT, dimension(0:2), parameter::coop_source_tau_n = (/ 1000, 850, 800 /)

  COOP_REAL, dimension(0:2), parameter::coop_source_k_weight = (/ 0.2d0, 0.15d0, 0.1d0 /)

  COOP_INT, dimension(0:2), parameter::coop_source_k_n = (/ 160, 120, 80 /)

  COOP_INT, dimension(0:2), parameter::coop_num_sources = (/ 4,  2,  2 /)

  COOP_REAL, parameter::coop_source_k_index = 0.55d0


  type coop_cosmology_firstorder_source
     COOP_INT::m = 0
     COOP_INT::ntau = 0
     COOP_INT::nk = 0
     COOP_INT::index_tc_max = 1
     COOP_INT, dimension(:),allocatable::index_tc_off
     COOP_INT, dimension(coop_pert_default_nq)::index_massivenu_on
     COOP_REAL, dimension(:),allocatable::k, tau, a, tauc, dk, dtau, chi  !!tau is conformal time, chi is comoving distance; in flat case, chi + tau = tau_0
     COOP_REAL, dimension(:,:,:),allocatable::s
   contains
     procedure::free => coop_cosmology_firstorder_source_free
  end type coop_cosmology_firstorder_source
  
  type, extends(coop_cosmology_background) :: coop_cosmology_firstorder
     logical::do_reionization = .true.
     COOP_REAL::zrecomb, distlss, tau0
     COOP_REAL::optre = 0.07d0
     COOP_REAL::zre = 8.d0
     COOP_REAL::deltaz = 0.5d0
     COOP_REAL::kMpc_pivot = 0.05d0
     COOP_INT::de_genre = COOP_PERT_NONE
     COOP_REAL::k_pivot
     COOP_REAL::dkappadtau_coef, ReionFrac, Omega_m, Omega_r, Omega_b, Omega_c, Omega_nu, Omega_g, a_eq, tau_eq, mnu_by_Tnu, Rbya

     type(coop_function)::Ps, Pt, Xe, ekappa, vis, Tb
     type(coop_cosmology_firstorder_source),dimension(0:2)::source
     COOP_INT::index_baryon, index_cdm, index_radiation, index_nu, index_massiveNu, index_de
     COOP_REAL, dimension(0:coop_pert_default_lmax, 0:coop_pert_default_mmax, 0:coop_pert_default_smax)::klms, klms_by_2lm1, klms_by_2lp1
     logical::klms_done = .false.
   contains
     procedure:: set_klms => coop_cosmology_firstorder_set_klms
     procedure:: set_source_tau => coop_cosmology_firstorder_set_source_tau
     procedure:: set_source_k => coop_cosmology_firstorder_set_source_k
     procedure:: set_power => coop_cosmology_firstorder_set_power
     procedure:: set_xe => coop_cosmology_firstorder_set_xe
     procedure:: set_zre_from_optre => coop_cosmology_firstorder_set_zre_from_optre
     procedure:: set_optre_from_zre => coop_cosmology_firstorder_set_optre_from_zre
     procedure::set_initial_conditions => coop_cosmology_firstorder_set_initial_conditions
     procedure:: xeofa => coop_cosmology_firstorder_xeofa
     procedure:: cs2bofa => coop_cosmology_firstorder_cs2bofa
     procedure:: Tbofa => coop_cosmology_firstorder_Tbofa
     procedure:: dxeda => coop_cosmology_firstorder_dxeda
     procedure:: dlnxedlna => coop_cosmology_firstorder_dlnxedlna
     procedure:: dkappada => coop_cosmology_firstorder_dkappada
     procedure:: dkappadtau => coop_cosmology_firstorder_dkappadtau
     procedure:: taucofa => coop_cosmology_firstorder_taucofa
     procedure:: dot_tauc => coop_cosmology_firstorder_dot_tauc
     procedure:: visofa => coop_cosmology_firstorder_visofa
     procedure:: ekappaofa => coop_cosmology_firstorder_ekappaofa
     procedure:: psofk => coop_cosmology_firstorder_psofk
     procedure:: ptofk => coop_cosmology_firstorder_ptofk
     procedure:: free => coop_cosmology_firstorder_free
     procedure::pert2source => coop_cosmology_firstorder_pert2source
     procedure::init_source => coop_cosmology_firstorder_init_source
     procedure::compute_source =>  coop_cosmology_firstorder_compute_source
     procedure::compute_source_k =>  coop_cosmology_firstorder_compute_source_k
  end type coop_cosmology_firstorder


contains



  subroutine coop_cosmology_firstorder_equations(n, tau, y, yp, cosmology, pert)
    COOP_INT n
    class(coop_cosmology_firstorder)::cosmology
    type(coop_pert_object)::pert
    COOP_REAL tau, y(0:n-1), yp(0:n-1)
    COOP_INT i, l, iq
    COOP_REAL a, aH, R, tauc, taucdot, aniso, ksq, ktauc, ktaucdot,  slip, P, aniso_dot, aHdot
    COOP_REAL rhoa2_g, pa2_g, rhoa2_b,  cs2b, rhoa2_c, rhoa2_nu, pa2_nu, rhoa2_mnu, pa2_mnu, rhoa2_sum, pa2_sum, Fmnu0(coop_pert_default_nq), Fmnu2(coop_pert_default_nq), Fmnu2_prime(coop_pert_default_nq), wmnu(coop_pert_default_nq), qbye(coop_pert_default_nq), pa2dot_g, pa2dot_nu, pa2dot_mnu, sumwmnu
    !!My PHI = Psi in Hu & White = Psi in Ma et al;
    !!My PSI = - Phi in Hu & White = Phi in Ma et al;
    !!My multipoles  = 4 * multipoles in Hu & White = (2l + 1) * multipoles in Ma et al
    !!My neutrino multipoles = (2l + 1) / (d ln f_0/d ln q) * Psi_l(q) in Ma et al
    yp(0) = 0
    a = cosmology%aoftau(tau)
    rhoa2_g = O0_RADIATION(cosmology)%rhoa2(a)
    pa2_g = O0_RADIATION(cosmology)%wofa(a) * rhoa2_g
    rhoa2_b = O0_BARYON(cosmology)%rhoa2(a)
    cs2b = cosmology%cs2bofa(a)
    rhoa2_c = O0_CDM(cosmology)%rhoa2(a)
    rhoa2_nu = O0_NU(cosmology)%rhoa2(a)
    pa2_nu = rhoa2_nu *  O0_NU(cosmology)%wofa(a)
    if(cosmology%index_massivenu .ne. 0)then
       rhoa2_mnu = O0_MASSIVENU(cosmology)%rhoa2(a)
       pa2_mnu = rhoa2_mnu *  O0_MASSIVENU(cosmology)%wofa(a)
       qbye =  coop_pert_default_q/sqrt(coop_pert_default_q**2 + (cosmology%mnu_by_Tnu * a)**2)
    else
       rhoa2_mnu = 0.d0
       pa2_mnu = 0.d0
       qbye = 1.d0
    endif
    wmnu =qbye*coop_pert_default_q_kernel
    sumwmnu = sum(wmnu)
    rhoa2_sum = rhoa2_g + rhoa2_b + rhoa2_c + rhoa2_nu +rhoa2_mnu
    pa2_sum = pa2_g + pa2_nu + pa2_mnu
    aH = sqrt(rhoa2_sum/3.d0)
    aHdot = -(rhoa2_sum+3.d0*pa2_sum)/6.d0

    pa2dot_nu = O0_NU(cosmology)%dpa2da(a)*aH*a
    pa2dot_g = O0_RADIATION(cosmology)%dpa2da(a)*aH*a
    if(pa2_mnu .eq. 0.d0)then
       pa2dot_mnu = 0.d0 
    else
       pa2dot_mnu = O0_MASSIVENU(cosmology)%dpa2da(a)*aH*a
    endif
    R = 0.75d0 * rhoa2_b/rhoa2_g
    ksq = pert%k ** 2
    tauc = cosmology%taucofa(a)
    taucdot = cosmology%dot_tauc(a)
    ktauc = pert%k * tauc
    ktaucdot = pert%k * taucdot
    select case(pert%m)
    case(0)

       O1_PSI_PRIME = O1_PSIDOT

       O1_DELTA_C_PRIME = - O1_V_C * pert%k + 3.d0 * O1_PSIDOT
       O1_DELTA_B_PRIME = - O1_V_B * pert%k + 3.d0 * O1_PSIDOT
       O1_T_PRIME(0) = - O1_T(1) * pert%k/3.d0 + 4.d0 * O1_PSIDOT 
       O1_NU_PRIME(0) = - O1_NU(1)*pert%k/3.d0 + 4.d0 * O1_PSIDOT 

       if(pert%tight_coupling)then
          pert%T%F(2) = (8.d0/9.d0)*ktauc * O1_T(1)
       else
          pert%T%F(2) = O1_T(2)
       endif
       aniso = pa2_g * pert%T%F(2) + pa2_nu * O1_NU(2)

       if(cosmology%index_massivenu .ne. 0)then
          if(pert%massivenu_iq_used .gt.0)then
             do iq = 1, pert%massivenu_iq_used
                Fmnu2(iq) = O1_MASSIVENU(2, iq)
                Fmnu0(iq) = O1_MASSIVENU(0, iq)
                O1_MASSIVENU_PRIME(0, iq) = - O1_NU(1)*pert%k/3.d0 * qbye(iq) + 4.d0 * O1_PSIDOT 
             enddo
             do iq = pert%massivenu_iq_used + 1, coop_pert_default_nq
                Fmnu2(iq) = O1_NU(2)
                Fmnu0(iq) = O1_NU(0)
             enddo
             aniso = aniso + pa2_mnu * sum(Fmnu2*wmnu)/sumwmnu
          else
             aniso = aniso + pa2_mnu * O1_NU(2)
          endif
       endif

       aniso = 0.6d0/ksq * aniso
       O1_PHI = O1_PSI - aniso

       !!velocities
       O1_V_C_PRIME = - aH * O1_V_C + pert%k * O1_PHI
       O1_NU_PRIME(1) = (O1_NU(0) + 4.d0*O1_PHI - 0.4d0 * O1_NU(2))*pert%k
       do iq=1, pert%massivenu_iq_used
          O1_MASSIVENU_PRIME(1, iq) = (O1_MASSIVENU(0, iq)*qbye(iq) + 4.d0*O1_PHI/qbye(iq) - 0.4d0 * O1_MASSIVENU(2, iq)*qbye(iq)) * pert%k
       enddo
       if(pert%tight_coupling)then
          slip = ktauc &
               * (- aH * O1_V_B + (cs2b * O1_DELTA_B - O1_T(0)/4.d0 +  pert%T%F(2)/10.d0)*pert%k ) &
               / ( pert%k * (R+1.d0)/R + ktaucdot )
          
       else
          slip = O1_V_B - O1_T(1)/4.d0
       endif
       O1_V_B_PRIME = -aH * O1_V_B + pert%k * (O1_PHI + cs2b * O1_DELTA_B) - slip/(R*tauc)
       O1_T_PRIME(1) = (O1_T(0) + 4.d0*O1_PHI - 0.4d0*pert%T%F(2))*pert%k + 4.d0*slip/tauc


       !!higher moments
       !!massless neutrinos
       do l = 2, pert%nu%lmax - 1
          O1_NU_PRIME(l) =  pert%k * (cosmology%klms_by_2lm1(l, 0, 0) *   O1_NU( l-1 ) - cosmology%klms_by_2lp1(l+1, 0, 0) *  O1_NU( l+1 ) )
       enddo
       O1_NU_PRIME(pert%nu%lmax) = pert%k *(pert%nu%lmax+0.5d0)/(pert%nu%lmax-0.5d0)*  O1_NU(pert%nu%lmax-1) &
            -  (pert%nu%lmax+1)* O1_NU(pert%nu%lmax)/tau

       if(cosmology%index_massivenu .ne. 0)then
          !!massive neutrinos
          do iq = 1, pert%massivenu_iq_used
             do l = 2, pert%massivenu(iq)%lmax - 1
                O1_MASSIVENU_PRIME(l, iq) = pert%k * qbye(iq) * (cosmology%klms_by_2lm1(l, 0, 0) * O1_MASSIVENU(l-1, iq) - cosmology%klms_by_2lp1(l+1, 0, 0) * O1_MASSIVENU(l+1,iq))
             enddo
             O1_MASSIVENU_PRIME(pert%massivenu(iq)%lmax, iq) = pert%k*pert%k * qbye(iq) * (pert%massivenu(iq)%lmax+0.5d0)/(pert%massivenu(iq)%lmax-0.5d0) *  O1_MASSIVENU(pert%nu%lmax-1, iq) &
                  -  (pert%nu%lmax+1)* O1_MASSIVENU(pert%nu%lmax, iq)/tau          
          enddo
       endif
       
       if(pert%tight_coupling)then
          pert%E%F(2) = -coop_sqrt6/4.d0 * pert%T%F(2)
          aniso_dot = pa2_g * (8.d0/9.d0)*(ktauc*O1_T_PRIME(1)+ktaucdot*O1_T(1)) + pa2dot_g * pert%T%F(2) + pa2_nu * O1_NU_PRIME(2) + pa2dot_nu*O1_NU(2)
          print*, a, tau
          print*, tau*aH
          print*, O1_NU(1)/3.d0*pert%k*tau, O1_NU(2), O1_NU_PRIME(2)/2.d0/aH
       else
          !!T
          P = (O1_T(2) - coop_sqrt6 * O1_E(2))/10.d0
          O1_T_PRIME(2) =  pert%k * (cosmology%klms_by_2lm1(2, 0, 0)*O1_T(1) - cosmology%klms_by_2lp1(3, 0, 0)*O1_T(3))  - (O1_T(2) - P)/tauc
          do l = 3, pert%T%lmax -1
             O1_T_PRIME(l) = pert%k * (cosmology%klms_by_2lm1(l, 0, 0)*O1_T(l-1) -cosmology%klms_by_2lp1(l+1, 0, 0)*O1_T(l+1))  - O1_T(l)/tauc
          enddo
          O1_T_PRIME(pert%T%lmax) =  pert%k *((pert%T%lmax+0.5d0)/(pert%T%lmax-0.5d0))*  O1_T(pert%T%lmax-1) &
            -  ((pert%T%lmax+1)/tau + 1.d0/tauc) * O1_T(pert%T%lmax)
          !!E
          O1_E_PRIME(2) = pert%k * ( - cosmology%klms_by_2lp1(3, 0, 0)*O1_E(3))  - (O1_E(2) + coop_sqrt6 * P)/tauc
          do l = 3, pert%E%lmax - 1
             O1_E_PRIME(l) = pert%k * (cosmology%klms_by_2lm1(l, 0, 2)*O1_E(l-1) - cosmology%klms_by_2lp1(l+1, 0, 2)*O1_E(l+1))  - O1_E(l)/tauc
          enddo
          O1_E_PRIME(pert%E%lmax) =  pert%k *( (pert%E%lmax+0.5d0)/(pert%E%lmax-0.5d0)) *  O1_E(pert%E%lmax-1) &
               -  ((pert%E%lmax+1)/tau + 1.d0/tauc) * O1_E(pert%E%lmax)
          aniso_dot = pa2_g * O1_T_PRIME(2) + pa2dot_g * O1_T(2) + pa2_nu * O1_NU_PRIME(2) + pa2dot_nu*O1_NU(2)
       endif
       if(cosmology%index_massivenu .ne. 0)then
          if(pert%massivenu_iq_used .le. 0)then
             aniso_dot = aniso_dot + pa2dot_mnu * O1_NU(2) + pa2_mnu + O1_NU_PRIME(2)
          else
             do iq = 1, pert%massivenu_iq_used
                Fmnu2_prime(iq) = O1_MASSIVENU_PRIME(2, iq)
             enddo
             do iq = pert%massivenu_iq_used + 1, coop_pert_default_nq
                Fmnu2_prime(iq) = O1_NU_PRIME(2)
             enddo
             aniso = aniso + pa2dot_mnu * sum(Fmnu2*wmnu)/sumwmnu + pa2_mnu * sum(Fmnu2_prime*wmnu)/sumwmnu 
          endif
       endif
       aniso_dot = 0.6d0/ksq * aniso_dot
       O1_PHI_PRIME = O1_PSI_PRIME - aniso_dot
       O1_PSIDOT_PRIME = - aH*(O1_PHI_PRIME + 3.d0*O1_PSI_PRIME) - 2.d0*(aHdot + aH**2)*O1_PHI - ksq/3.d0*(O1_PSI+aniso) + (rhoa2_b * O1_DELTA_B * (cs2b - 1.d0/3.d0) + rhoa2_c*O1_DELTA_C*(-1.d0/3.d0))/2.d0

       if(cosmology%index_massivenu .ne. 0)then
          if(pert%massivenu_iq_used .gt. 0)then
             O1_PSIDOT_PRIME =  O1_PSIDOT_PRIME + (pa2_mnu*sum(Fmnu0*wmnu)/sumwmnu - rhoa2_mnu/3.d0  * sum(Fmnu0*coop_pert_default_q_kernel/qbye)/sum(coop_pert_default_q_kernel/qbye) )/2.d0
          endif
        endif
        select case(pert%de%genre)
        case(COOP_PERT_NONE)
           !!do nothing
        case default
           call coop_tbw("de perturbations not written")
        end select
    case(1)
       call coop_tbw("vector equations not written")
    case(2)
       call coop_tbw("tensor equations not written")
    case default
       call coop_return_error("firstorder_equations", "Unknown m = "//trim(coop_num2str(pert%m)), "stop")
    end select
  end subroutine coop_cosmology_firstorder_equations



  subroutine coop_cosmology_firstorder_pert2source(this, pert, source, itau, ik)
    class(coop_cosmology_firstorder)::this
    type(coop_pert_object)::pert
    COOP_INT itau, ik
    type(coop_cosmology_firstorder_source)::source
    select case(source%m)
    case(0)
       source%s(itau,  ik, 1) = pert%O1_PSI
       source%s(itau,  ik, 2) = pert%O1_PSIDOT
       source%s(itau,  ik, 3) = pert%O1_T(0)
    case(1)
       call coop_tbw("vector source to be done")
    case(2)
       call coop_tbw("tensor source to be done")
    end select
  end subroutine coop_cosmology_firstorder_pert2source


  subroutine coop_cosmology_firstorder_set_initial_conditions(this, pert, m, k, tau)
    class(coop_cosmology_firstorder)::this
    class(coop_pert_object)::pert
    COOP_INT :: m 
    COOP_REAL tau, k, Rnu
    pert%k = k
    pert%tight_coupling = .true.
    call pert%init(m = m, nu_mass = this%mnu_by_Tnu, de_genre = this%de_genre)
    call pert%set_zero()
    select case(trim(pert%initial_conditions))
    case("adiabatic")
       select case(pert%m)
       case(0)
          Rnu = this%Omega_nu / this%Omega_r
          pert%O1_Phi = coop_primordial_zeta_norm/(1.5d0 + 0.4d0*Rnu)
          pert%O1_Phidot = 0.d0
          pert%O1_PSI =  (1.+0.4*Rnu)*pert%O1_Phi
          pert%O1_PSIDOT = 0.d0

          pert%O1_DELTA_B = -1.5d0*pert%O1_Phi
          pert%O1_DELTA_C = pert%O1_DELTA_B
          pert%O1_T(0) =  -2.d0*pert%O1_Phi
          pert%O1_NU(0) = -2.d0*pert%O1_T(0)
          

          pert%O1_V_B = pert%O1_Phi/2.d0*k*tau
          pert%O1_V_C = pert%O1_V_B
          pert%O1_NU(1) = 4.d0*pert%O1_V_C
          pert%O1_T(1) = pert%O1_NU(1)

          pert%O1_NU(2) = (2.d0/3.d0) * pert%O1_Phi * (k*tau)**2
       case(1)
          call coop_tbw("vector initialization")
       case(2)
          call coop_tbw("tensor initialization")
       end select
    case default
       stop "unknown initial conditions"
    end select
  end subroutine coop_cosmology_firstorder_set_initial_conditions




  subroutine coop_cosmology_firstorder_set_klms(this)
    class(coop_cosmology_firstorder)::this
    COOP_INT k, l, m, s
    if(this%klms_done) return
    do s = 0, coop_pert_default_smax
       do m = 0,coop_pert_default_mmax
          do l = 0, coop_pert_default_lmax
             if(m.ge.l .or. s.ge.l)then
                this%klms(l,m,s) = 0.d0
             else
                this%klms(l, m, s) = sqrt(dble(l**2-m**2)*dble(l**2-s**2))/l
             endif
             this%klms_by_2lp1(l, m, s) = this%klms(l, m, s)/(2*l+1.d0)
             this%klms_by_2lm1(l, m, s) = this%klms(l, m, s)/(2*l-1.d0)
          enddo
       enddo
    enddo
    this%klms_done = .true.
  end subroutine coop_cosmology_firstorder_set_klms


  function coop_cosmology_firstorder_cs2bofa(this, a) result(cs2bofa)
    class(coop_cosmology_firstorder)::this
    COOP_REAL a, cs2bofa
    cs2bofa = this%species(this%index_baryon)%fcs2%eval(a)
  end function coop_cosmology_firstorder_cs2bofa

  function coop_cosmology_firstorder_Tbofa(this, a) result(Tbofa)
    class(coop_cosmology_firstorder)::this
    COOP_REAL a, Tbofa
    Tbofa = this%Tb%eval(a)
  end function coop_cosmology_firstorder_Tbofa


  subroutine coop_cosmology_firstorder_source_free(this)
    class(coop_cosmology_firstorder_source)::this
    if(allocated(this%k))deallocate(this%k, this%dk, this%index_tc_off)
    this%nk = 0
    if(allocated(this%tau))deallocate(this%tau, this%chi, this%dtau, this%a, this%tauc)
    this%ntau = 0
    if(allocated(this%s))deallocate(this%s)
  end subroutine coop_cosmology_firstorder_source_free

  subroutine coop_cosmology_firstorder_free(this)
    class(coop_cosmology_firstorder)::this
    COOP_INT m, i
    do m=0,2
       call this%source(m)%free()
    enddo
    call this%Ps%free()
    call this%Pt%free()
    call this%Xe%free()
    call this%eKappa%free()
    call this%vis%free()
    call this%Tb%free()
    call this%fdis%free
    call this%ftime%free
    call this%faoftau%free
    do i= 1, this%num_species
       call this%species(i)%free
    enddo
    this%num_species = 0
    this%omega_k_value = 1.
    this%need_setup_background = .true.
  end subroutine coop_cosmology_firstorder_free


  subroutine coop_cosmology_firstorder_set_power(this, power, args)
!!the format of primordial power
!!subroutine power(k/k_pivot, ps, pt, cosmology, args)
!!input kMpc and args, output ps and pt    
    integer,parameter::n = 4096
    class(coop_cosmology_firstorder)::this
    external power
    type(coop_arguments) args
    COOP_REAL k(n), ps(n), pt(n)
    COOP_INT i
    this%k_pivot = this%kMpc_pivot/this%H0Mpc()
    call coop_set_uniform(n, k, coop_power_lnk_min, coop_power_lnk_max)
    k = exp(k)
    do i=1, n
       call power(k(i)/this%k_pivot, ps(i), pt(i), this, args)
    end do
    call this%ps%init(n = n, xmin = k(1), xmax = k(n), f = ps, xlog = .true., ylog = .true., fleft = ps(1), fright = ps(n), check_boundary = .false.)
    call this%pt%init(n = n, xmin = k(1), xmax = k(n), f = pt, xlog = .true., ylog = .true., fleft = pt(1), fright = pt(n), check_boundary = .false.)

  end subroutine coop_cosmology_firstorder_set_power

  subroutine coop_cosmology_firstorder_set_xe(this)
    class(coop_cosmology_firstorder)::this
    integer i, m
    integer, parameter::n = 4096
    COOP_REAL ekappa(n), a(n), dkappadinva(n), dkappadtau(n), vis(n), acal
    call this%set_klms()
    this%index_baryon = this%index_of("Baryon", .true.)
    this%index_radiation = this%index_of("Radiation", .true.)
    this%index_cdm = this%index_of("CDM", .true.)
    this%index_Nu = this%index_of("Massless Neutrinos", .true.)
    this%index_massiveNu = this%index_of("Massive Neutrinos")
    this%index_de = this%index_of("Dark Energy", .true.)
    
    acal  = coop_visibility_amin

    !!asymptotic Omega_b and Omega_c
    this%Omega_b = O0_BARYON(this)%Omega *  O0_BARYON(this)%rhoa3_ratio(acal)
    this%Omega_c = O0_CDM(this)%Omega * O0_CDM(this)%rhoa3_ratio(acal)

    this%Omega_m = this%Omega_b + this%Omega_c
    this%Omega_g = O0_RADIATION(this)%Omega * O0_RADIATION(this)%rhoa4_ratio(acal)
    if(this%index_massiveNu .ne. 0)then
       call coop_fermion_get_lnam(log(O0_MASSIVENU(this)%Omega/O0_MASSIVENU(this)%Omega_massless), this%mnu_by_Tnu)
       this%mnu_by_Tnu = exp(this%mnu_by_Tnu)
       this%Omega_nu = O0_NU(this)%Omega * O0_NU(this)%rhoa4_ratio(acal) + O0_MASSIVENU(this)%Omega * O0_MASSIVENU(this)%rhoa4_ratio(acal)
    else
       this%Omega_nu = O0_NU(this)%Omega * O0_NU(this)%rhoa4_ratio(acal)
       this%mnu_by_Tnu = 0.d0
    endif
    this%Rbya = (3.d0/4.d0)*this%Omega_b/this%Omega_g

    this%Omega_r = this%Omega_g + this%Omega_nu 

    this%a_eq = this%Omega_r / this%Omega_m
    this%tau_eq = this%conformal_time(this%a_eq)

    this%dkappadtau_coef = this%species(this%index_baryon)%Omega * this%h() * coop_SI_sigma_thomson * (coop_SI_rhocritbyh2/coop_SI_c**2) * coop_SI_hbyH0 * coop_SI_c/ coop_SI_m_H * (1.d0 - this%YHe())
    if(this%do_reionization)then
       this%reionFrac = 1.d0 + this%YHe()/(coop_m_He_by_m_H * (1.d0-  this%YHe()))
    else
       this%reionFrac = 0.d0
    endif
    call this%set_zre_from_optre()

    call coop_recfast_get_xe(this, this%xe, this%Tb, this%reionFrac, this%zre, this%deltaz)
    call coop_set_uniform(n, a, log(coop_visibility_amin), log(coop_scale_factor_today))
    a = exp(a)
    !$omp parallel do
    do i=1, n
       dkappadtau(i) = this%dkappadtau(a(i))
       dkappadinva(i)= - dkappadtau(i)/this%Hratio(a(i))
    enddo
    !$omp end parallel do
    ekappa(n) = 0.d0
    do i=n-1,1, -1
       ekappa(i) = ekappa(i+1) + (dkappadinva(i+1) + dkappadinva(i))*(0.5d0/a(i) - 0.5d0/a(i+1))
    enddo
    ekappa = exp(ekappa)
    vis = dkappadtau*ekappa
    call this%ekappa%init(n = n, xmin = a(1), xmax = a(n), f = ekappa, xlog = .true., ylog = .true., fleft = ekappa(1), fright = ekappa(n), check_boundary = .false.)
    call this%vis%init(n = n, xmin = a(1), xmax = a(n), f = vis, xlog = .true., ylog = .true., fleft = vis(1), fright = vis(n), check_boundary = .false.)
    this%zrecomb = 1.d0/this%vis%maxloc() - 1.d0
    this%distlss = this%comoving_distance(1.d0/(1.d0+this%zrecomb))
    this%tau0 = this%conformal_time(coop_scale_factor_today)
    
  end subroutine coop_cosmology_firstorder_set_xe

  function coop_cosmology_firstorder_xeofa(this, a) result(xe)
    COOP_REAL a, xe
    class(coop_cosmology_firstorder)::this
    if(a .lt. coop_recfast_a_initial)then
       xe = this%xe%fleft
    else
       xe = this%xe%eval(a)
    endif
  end function coop_cosmology_firstorder_xeofa

  function coop_cosmology_firstorder_dxeda(this, a) result(dxeda)
    COOP_REAL a, dxeda
    class(coop_cosmology_firstorder)::this
    if(a .lt. coop_recfast_a_initial)then
       dxeda = 0.d0 
    else
       dxeda = this%xe%derivative(a)
    endif
  end function coop_cosmology_firstorder_dxeda

  function coop_cosmology_firstorder_dlnxedlna(this, a) result(dlnxedlna)
    COOP_REAL a, dlnxedlna
    class(coop_cosmology_firstorder)::this
    if(a .lt. coop_recfast_a_initial)then
       dlnxedlna = 0
    else
       dlnxedlna = this%xe%derivative_bare(log(a))
    endif
  end function coop_cosmology_firstorder_dlnxedlna



  function coop_cosmology_firstorder_ekappaofa(this, a) result(ekappa)
    COOP_REAL a, ekappa
    class(coop_cosmology_firstorder)::this
    if(a .lt. 1.1d-4)then
       ekappa = 0.d0
    else
       ekappa = this%ekappa%eval(a)
    endif    
  end function coop_cosmology_firstorder_ekappaofa

  function coop_cosmology_firstorder_visofa(this, a) result(vis)
    COOP_REAL a, vis
    class(coop_cosmology_firstorder)::this
    if(a .lt. 1.1d-4)then
       vis = 0.d0
    else
       vis = this%vis%eval(a)
    endif    
  end function coop_cosmology_firstorder_visofa


  function coop_cosmology_firstorder_dkappadtau(this, a) result(dkappadtau)
    class(coop_cosmology_firstorder)::this    
    COOP_REAL a, dkappadtau
    dkappadtau =  this%dkappadtau_coef * this%xeofa(a)/a**2
  end function coop_cosmology_firstorder_dkappadtau


  function coop_cosmology_firstorder_taucofa(this, a) result(tauc)
    class(coop_cosmology_firstorder)::this    
    COOP_REAL a, tauc
    tauc =   (this%dkappadtau_coef * this%xeofa(a))
    if(tauc .gt. 0.d0)then
       tauc =   a**2/tauc
    else
       tauc = 1.d99
    endif
  end function coop_cosmology_firstorder_taucofa


  function coop_cosmology_firstorder_dot_tauc(this, a) result(dtauc)
    class(coop_cosmology_firstorder)::this    
    COOP_REAL a, dtauc
    dtauc =  a*(2.d0 - this%dlnxedlna(a))/this%xeofa(a)/this%dkappadtau_coef*this%Hasq(a)
  end function coop_cosmology_firstorder_dot_tauc



  function coop_cosmology_firstorder_dkappada(this, a) result(dkappada)
    class(coop_cosmology_firstorder)::this    
    COOP_REAL a, dkappada
    dkappada = this%dkappadtau(a) / this%dadtau(a)
  end function coop_cosmology_firstorder_dkappada


  subroutine coop_cosmology_firstorder_set_source_tau(this, source, n, weight)
    COOP_REAL, parameter::incr = 1.05d0
    COOP_INT::nbuffer = 10
    COOP_INT n, i, j
    COOP_REAL  top(n), dtop, amax, amin, amid, topmax, topmin, topmid, weight, taucut, acut
    class(coop_cosmology_firstorder_source)::source
    class(coop_cosmology_firstorder)::this
    if(n.lt. 30) stop "n is too small for set_source_tau"
    source%ntau = n
    if(allocated(source%a))then
       if(size(source%a).ne.n)then
          deallocate(source%a, source%tau, source%chi, source%dtau, source%tauc)
          allocate(source%a(n), source%tau(n), source%dtau(n), source%chi(n), source%tauc(n))
       endif
    else
       allocate(source%a(n), source%tau(n), source%dtau(n), source%chi(n),  source%tauc(n))
    endif
    dtop = this%tau0 / (n - nbuffer)
    do i= nbuffer+1, n
       top(i) = dtop*(i-nbuffer-0.5d0)
    enddo
    do i=nbuffer, 1, -1
       top(i) = top(i+1)*0.001d0
    enddo
    amin = coop_visibility_amin
    topmin = this%ekappaofa(amin)**weight*this%conformal_time(amin)
    do i=1, n
       amax = amin*incr
       topmax = this%ekappaofa(amax)**weight*this%conformal_time(amax)
       do while( topmax .lt. top(i))
          amin = amax
          topmin = topmax
          amax = amax*incr
          if(amax .ge. coop_scale_factor_today)then
             amax = coop_scale_factor_today
             topmax = this%tau0
             exit
          endif
          topmax = this%ekappaofa(amax)**weight*this%conformal_time(amax)       
       enddo
       do while((amax-amin)/amin .gt. 1.d-7)
          amid = (amin+amax)/2.d0
          topmid = this%ekappaofa(amid)**weight*this%conformal_time(amid)       
          if(topmid .gt. top(i))then
             amax = amid
             topmax = topmid
          else
             amin = amid
             topmin = topmid
          endif
       enddo
       source%a(i) = (amin+amax)/2.d0
       source%tau(i) = this%conformal_time(source%a(i))
       if(i.gt.nbuffer+1)source%dtau(i) = dtop/((weight*this%dkappadtau(source%a(i))*source%tau(i)+1.d0)*this%ekappaofa(source%a(i))**weight)
       source%tauc(i) = this%taucofa(source%a(i))
    enddo
    source%dtau(1) = (source%tau(1)+source%tau(2))/2.d0
    do i=2, nbuffer+1
       source%dtau(i) = (source%tau(i+1) - source%tau(i-1))/2.d0
    enddo
    source%chi = this%tau0 - source%tau
    source%index_tc_max = 1
    do while( source%index_tc_max .lt. source%ntau - 1 .and. source%tauc(source%index_tc_max+1) .le. 0.03 * source%tau(source%index_tc_max+1) )
       source%index_tc_max = source%index_tc_max+1
    enddo
  end subroutine coop_cosmology_firstorder_set_source_tau
  
  subroutine coop_cosmology_firstorder_set_source_k(this, source, n, weight)
    class(coop_cosmology_firstorder)::this
    class(coop_cosmology_firstorder_source)::source
    COOP_INT n, i, iq
    COOP_REAL weight, kop(n), dkop, kopmin, kopmax, kmin, kmax, kalpha
    this%k_pivot = this%kMpc_pivot/this%H0Mpc()
    source%nk = n
    if(allocated(source%k))then
       if(size(source%k).ne.n)then
          deallocate(source%k, source%dk, source%index_tc_off)
          allocate(source%k(n), source%dk(n), source%index_tc_off(n))
       endif
    else
       allocate(source%k(n), source%dk(n), source%index_tc_off(n))
    endif
    kmin = 0.3d0/this%distlss
    select case(source%m)
    case(0)
       kmax = 6.2d3/this%distlss
    case(1)
       kmax = 4.2d3/this%distlss
    case(2)
       kmax = 2.2d3/this%distlss
    case default
       write(*,*) "Error: m = ", source%m
       stop "source spin must be 0, 1, 2"
    end select

    kopmin = kmin**coop_source_k_index*weight + log(kmin)
    kopmax = kmax**coop_source_k_index*weight + log(kmax)
    dkop = (kopmax-kopmin)/(n-1)
    call coop_set_uniform(n, kop, kopmin, kopmax)
    !$omp parallel do private(kalpha)
    do i=1, n
       call coop_source_kop2k_noindex(kop(i)*coop_source_k_index, weight*coop_source_k_index, kalpha)
       source%k(i) = kalpha**(1.d0/coop_source_k_index)
       source%dk(i) = dkop * source%k(i) / (1.d0 + weight * coop_source_k_index * kalpha )
    enddo
    !$omp end parallel do

    !!set index_tc_off
    if(source%ntau .gt. 0 .and. allocated(source%tauc))then
       source%index_tc_off = 1
       do i = 1, n
          do while(source%index_tc_off(i) .lt. source%index_tc_max .and. source%tauc(source%index_tc_off(i)+1)*source%k(i) .le. coop_cosmology_firstorder_tc_cutoff)
             source%index_tc_off(i)= source%index_tc_off(i)+1
          enddo
       enddo
    else
       call coop_return_error("set_source_k", "You need to call set_source_tau before calling set_source_k", "stop")
    endif
    
    !!set index_massivenu_on
    if(this%mnu_by_Tnu*source%a(source%ntau) .lt. coop_pert_massivenu_threshold(1))then
       source%index_massivenu_on = source%ntau + 1
       return
    endif
    source%index_massivenu_on(coop_pert_default_nq) = source%ntau + 1
    iq = coop_pert_default_nq
    do while(this%mnu_by_Tnu*source%a(source%index_massivenu_on(iq)-1) .gt. coop_pert_massivenu_threshold(iq))
       source%index_massivenu_on(iq) = source%index_massivenu_on(iq)-1
       if(source%index_massivenu_on(iq) .le. 1)exit
    enddo
    do iq = coop_pert_default_nq - 1, 1, -1
       source%index_massivenu_on(iq) =  source%index_massivenu_on(iq+1)
       do while(this%mnu_by_Tnu*source%a(source%index_massivenu_on(iq)-1) .gt. coop_pert_massivenu_threshold(iq))
          source%index_massivenu_on(iq) = source%index_massivenu_on(iq)-1
          if(source%index_massivenu_on(iq) .le. 1)exit
       enddo
    enddo
    
  end subroutine coop_cosmology_firstorder_set_source_k


  subroutine coop_source_kop2k_noindex(kop, weight, k)
    !!solve the equation k * weight + log(k) = kop
    COOP_REAL kop, weight, k, kmin, kmax, kmid
    if(kop .gt. weight)then
       kmax = kop/weight
       kmin = (kop+1.d0)/(weight+1.d0)
    else
       kmax = exp(kop)
       kmin = max(kop/weight, 1.d-3)
    endif
    do while((kmax-kmin)/kmin .gt. 1.d-8)
       kmid = (kmax+kmin)/2.d0
       if(kmid*weight + log(kmid) .gt. kop)then
          kmax = kmid
       else
          kmin = kmid
       endif
    enddo
    k = (kmax+kmin)/2.d0
  end subroutine coop_source_kop2k_noindex
  


  function coop_cosmology_firstorder_psofk(this, k) result(ps)
    COOP_REAL k, ps
    class(coop_cosmology_firstorder)::this
    ps = this%ps%eval(k)
  end function coop_cosmology_firstorder_psofk


  function coop_cosmology_firstorder_ptofk(this, k) result(pt)
    COOP_REAL k, pt
    class(coop_cosmology_firstorder)::this
    pt = this%ps%eval(k)
  end function coop_cosmology_firstorder_ptofk


  subroutine coop_cosmology_firstorder_set_zre_from_optre(this)
    class(coop_cosmology_firstorder)::this
    COOP_REAL zremin, zremax, zremid
    COOP_REAL optremin, optremax, optremid, optre_wanted
    integer iloop
    if(this%ReionFrac .le. 0.d0)then
       this%zre = -1.d0
    endif
    optre_wanted = this%optre

    zremin = 0.d0
    this%zre = zremin
    call this%set_optre_from_zre()
    optremin = this%optre

    zremax = 20.d0
    this%zre = zremax
    call this%set_optre_from_zre()    
    optremax = this%optre

    do while(optre_wanted .gt. optremax)
       zremin = zremax
       optremin = optremax
       zremax = zremax * 1.3d0
       this%zre = zremax
       call this%set_optre_from_zre()    
       optremax = this%optre
    enddo
    iloop = 0
    do while((optremax - optremin) .gt. 1.d-5 .and. iloop .lt. 25)
       zremid = (zremin + zremax)/2.d0
       this%zre = zremid
       call this%set_optre_from_zre()
       optremid = this%optre
       if(optremid .gt. optre_wanted)then
          zremax = zremid
          optremax = optremid
       else
          zremin = zremid
          optremin = optremid
       endif
       iloop = iloop + 1
    enddo
    this%optre = optre_wanted
  end subroutine coop_cosmology_firstorder_set_zre_from_optre

  subroutine coop_cosmology_firstorder_set_optre_from_zre(this)
    class(coop_cosmology_firstorder)::this
    COOP_REAL :: optre, ReionFrac, zre, delta
    type(coop_arguments):: args
    if(this%ReionFrac .le. 0.d0)then
       this%optre = 0.d0
    else
       this%optre = coop_integrate(coop_optre_int, 0.d0, this%zre + this%deltaz * 10.d0, this, args, COOP_REAL_OF(1.e-7))
    endif
  end subroutine coop_cosmology_firstorder_set_optre_from_zre

  function coop_optre_int(z, cosmology, args) result(dkappadz)
    COOP_REAL z, dkappadz
    class(coop_cosmology_firstorder) cosmology
    type(coop_arguments)::args
    dkappadz = cosmology%dkappadtau_coef * coop_reionization_xe(z, cosmology%reionFrac, cosmology%zre, cosmology%deltaz) / cosmology%Hasq(1.d0/(1.d0+z)) 
  end function coop_optre_int





  subroutine coop_cosmology_firstorder_init_source(this, m)
    class(coop_cosmology_firstorder)::this
    COOP_INT :: m
    this%source(m)%m = m
    call this%set_source_tau(this%source(m), coop_source_tau_n(m), coop_source_tau_weight(m))
    call this%set_source_k(this%source(m), coop_source_k_n(m), coop_source_k_weight(m))
    if(allocated(this%source(m)%s))then
       if(size(this%source(m)%s, 1) .ne. this%source(m)%ntau .or. size(this%source(m)%s, 2) .ne. this%source(m)%nk)then
          deallocate(this%source(m)%s)
          allocate(this%source(m)%s(this%source(m)%ntau, this%source(m)%nk, coop_num_sources(m)))
       endif
    else
       allocate(this%source(m)%s(this%source(m)%ntau, this%source(m)%nk, coop_num_sources(m)))
    endif
  end subroutine coop_cosmology_firstorder_init_source


  subroutine coop_cosmology_firstorder_compute_source(this, m)
    class(coop_cosmology_firstorder)::this
    COOP_INT m, ik
    call this%init_source(m)
   ! !$omp parallel do
    do ik = 1, this%source(m)%nk
       call this%compute_source_k(this%source(m), ik)
    enddo
  ! !$omp end parallel do
  end subroutine coop_cosmology_firstorder_compute_source


  subroutine coop_cosmology_firstorder_compute_source_k(this, source, ik)
    class(coop_cosmology_firstorder)::this
    type(coop_cosmology_firstorder_source)::source
    COOP_INT ik, nvars, nw, itau, iq
    type(coop_pert_object) pert
    COOP_REAL, dimension(:,:),allocatable::w
    COOP_REAL c(24)
    COOP_INT ind, i
    COOP_REAL tau_ini, tau
    tau_ini = min(coop_initial_condition_epsilon/source%k(ik), this%conformal_time(this%a_eq*coop_initial_condition_epsilon))
    call this%set_initial_conditions(pert, m = source%m, k = source%k(ik), tau = tau_ini)
    call coop_cosmology_firstorder_equations(pert%ny+1, tau, pert%y, pert%yp, this, pert)
    do i=0, pert%ny
       print*, pert%y(i), pert%yp(i)
    enddo
    stop
    ind = 1
    c = 0.d0
    nvars = pert%ny + 1
    nw = nvars
    allocate(w(nw, 9))
    w = 0.d0
    pert%tau = tau_ini
    iq = 1
    do itau = 1, source%ntau
       call coop_dverk_firstorder(nvars, coop_cosmology_firstorder_equations, this, pert, pert%tau,   pert%y(0:pert%ny), source%tau(itau),  coop_cosmology_firstorder_ode_accuracy, ind, c, nw, w)
       call this%pert2source(pert, source, itau, ik)
       if(itau .eq. source%index_tc_off(ik))then
          call pert%save_ode()
          if(pert%m .eq. 0)then  !!set tight coupling approximations
             pert%T%F(2) = (8.d0/9.d0) * pert%k * this%taucofa(source%a(itau)) * pert%T%F(1)
             pert%E%F(2) = -(coop_sqrt6/4.d0)*pert%T%F(2)
          endif
          pert%tight_coupling = .false.
          call pert%init(m = source%m, nu_mass = this%mnu_by_Tnu, de_genre = this%de_genre, a = source%a(itau))
          call pert%restore_ode()
          ind = 1
          c = 0.d0
          deallocate(w)
          nvars = pert%ny + 1
          nw = nvars
          allocate(w(nw, 9))
       endif
       if(iq.le.coop_pert_default_nq)then
          do while(itau .eq. source%index_massivenu_on(iq))
             call pert%save_ode()
             call pert%init(m = source%m, nu_mass = this%mnu_by_Tnu, de_genre = this%de_genre, a = source%a(itau))             
             call pert%restore_ode()
             ind = 1
             c = 0.d0
             deallocate(w)
             nvars = pert%ny + 1
             nw = nvars
             allocate(w(nw, 9))
             iq = iq + 1
             if(iq .gt.coop_pert_default_nq) exit
          enddo
       endif
    enddo
    deallocate(w)
    call pert%free()
  end subroutine coop_cosmology_firstorder_compute_source_k



  subroutine coop_dverk_firstorder(n, fcn, cosmology, pert, x, y, xend, tol, ind, c, nw, w)
    implicit COOP_REAL (a-h,o-z)
    implicit COOP_INT (i-n)
    class(coop_cosmology_firstorder) cosmology
    type(coop_pert_object) pert
#define DVERK_ARGUMENTS ,cosmology,pert
#include "dverk.h"    
#undef DVERK_ARGUMENTS
  end subroutine coop_dverk_firstorder



end module coop_firstorder_mod
