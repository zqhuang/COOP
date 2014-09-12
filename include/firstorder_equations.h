  subroutine coop_cosmology_firstorder_equations(n, lna, y, yp, cosmology, pert)
    COOP_INT n
    type(coop_cosmology_firstorder)::cosmology
    type(coop_pert_object)::pert
    COOP_REAL lna, y(0:n-1), yp(0:n-1)
    COOP_INT i, l, iq
    COOP_REAL a, aniso, kbyaHsq, ktauc, ktaucdot, P, aniso_prime, kbyaH, aHtauc, aHtau, ksq, aHsq, uterm, vterm
    COOP_REAL :: pa2pr_g, pa2pr_nu
    COOP_REAL, dimension(coop_pert_default_nq)::Fmnu2_prime, Fmnu2, Fmnu0, qbye, wrho, wp, wrho_prime, wp_prime
    !!My PHI = Psi in Hu & White = Psi in Ma et al;
    !!My PSI = - Phi in Hu & White = Phi in Ma et al;
    !!My multipoles  = 4 * multipoles in Hu & White = (2l + 1) * multipoles in Ma et al
    !!My neutrino multipoles = (2l + 1) / (d ln f_0/d ln q) * Psi_l(q) in Ma et al
    !!finally, my time variable is log(a)
    yp(0) = 0
    ksq = pert%k**2
    a = exp(lna)
    pert%a = a
    pert%rhoa2_g = O0_RADIATION(cosmology)%rhoa2(a)
    pert%pa2_g = O0_RADIATION(cosmology)%wofa(a) * pert%rhoa2_g
    pert%rhoa2_b = O0_BARYON(cosmology)%rhoa2(a)
    pert%cs2b = cosmology%cs2bofa(a)
    pert%rhoa2_c = O0_CDM(cosmology)%rhoa2(a)
    pert%rhoa2_nu = O0_NU(cosmology)%rhoa2(a)
    pert%pa2_nu = pert%rhoa2_nu *  O0_NU(cosmology)%wofa(a)
    pert%rhoa2_de = O0_DE(cosmology)%rhoa2(a)
    pert%pa2_de = O0_DE(cosmology)%wofa(a)* pert%rhoa2_de
    if(cosmology%index_massivenu .ne. 0)then
       pert%rhoa2_mnu = O0_MASSIVENU(cosmology)%rhoa2(a)
       pert%pa2_mnu = pert%rhoa2_mnu *  O0_MASSIVENU(cosmology)%wofa(a)
       qbye =  coop_pert_default_q/sqrt(coop_pert_default_q**2 + (cosmology%mnu_by_Tnu * a)**2)
       wrho = coop_pert_default_q_kernel/qbye *  pert%num_mnu_ratio
       wrho_prime =   (cosmology%mnu_by_Tnu * a)**2/(coop_pert_default_q**2 + (cosmology%mnu_by_Tnu * a)**2) * wrho

       wp = qbye * coop_pert_default_q_kernel * pert%num_mnu_ratio
       wp_prime = -(cosmology%mnu_by_Tnu * a)**2/(coop_pert_default_q**2 + (cosmology%mnu_by_Tnu * a)**2) * wp

    else
       pert%rhoa2_mnu = 0.d0
       pert%pa2_mnu = 0.d0
    endif

    pert%rhoa2_sum = pert%rhoa2_g + pert%rhoa2_b + pert%rhoa2_c + pert%rhoa2_nu +pert%rhoa2_mnu+pert%rhoa2_de 
    pert%pa2_sum = pert%pa2_g + pert%pa2_nu + pert%pa2_mnu + pert%pa2_de
    aHsq = (pert%rhoa2_sum + cosmology%Omega_k())/3.d0

    pert%aH = sqrt(aHsq)

    pert%daHdtau = -(pert%rhoa2_sum+3.d0*pert%pa2_sum)/6.d0

    pa2pr_nu = O0_NU(cosmology)%dpa2da(a)*a
    pa2pr_g = O0_RADIATION(cosmology)%dpa2da(a)*a

    pert%R = 0.75d0 * pert%rhoa2_b/pert%rhoa2_g
    pert%tauc = cosmology%taucofa(a)
    pert%taucdot = cosmology%dot_tauc(a)

    ktauc = pert%k * pert%tauc
    ktaucdot = pert%k * pert%taucdot
    kbyaH  = pert%k/pert%aH
    kbyaHsq = kbyaH**2
    aHtauc = pert%aH * pert%tauc

    pert%tau = cosmology%tauofa(a)
    aHtau = pert%aH*pert%tau
 
    select case(pert%m)
    case(0)

       O1_PSI_PRIME = O1_PSIPR

       O1_DELTA_C_PRIME = - O1_V_C * kbyaH + 3.d0 * O1_PSI_PRIME
       O1_DELTA_B_PRIME = - O1_V_B * kbyaH + 3.d0 * O1_PSI_PRIME
       O1_T_PRIME(0) = - O1_T(1) * kbyaH/3.d0 + 4.d0 * O1_PSI_PRIME 
       O1_NU_PRIME(0) = - O1_NU(1) * kbyaH/3.d0 + 4.d0 * O1_PSI_PRIME

       if(pert%tight_coupling)then
          pert%T%F(2) = (8.d0/9.d0)*ktauc * O1_T(1)
          Uterm = - pert%aH * O1_V_B + (pert%cs2b * O1_DELTA_B - O1_T(0)/4.d0 +  pert%T%F(2)/10.d0)*pert%k
          vterm = pert%k * (pert%R+1.d0)/pert%R + ktaucdot 
          pert%slip = ktauc &
               * (Uterm - ktauc*(Uterm*(-pert%aH*pert%k/pert%R + pert%daHdtau/pert%aH*ktaucdot)/vterm)/vterm)/vterm   !!v_b - v_g accurate to (k tau_c)^2 
       else
          pert%T%F(2) = O1_T(2)
          pert%slip = O1_V_B - O1_T(1)/4.d0
       endif
       aniso = pert%pa2_g * pert%T%F(2) + pert%pa2_nu * O1_NU(2)

       if(cosmology%index_massivenu .ne. 0)then
          do iq = 1, pert%massivenu_iq_used
             Fmnu2(iq) = O1_MASSIVENU(2, iq)
             Fmnu0(iq) = O1_MASSIVENU(0, iq)
             O1_MASSIVENU_PRIME(0, iq) = - O1_NU(1)*kbyaH/3.d0 * qbye(iq) + 4.d0 * O1_PSI_PRIME 
          enddo
          do iq = pert%massivenu_iq_used + 1, coop_pert_default_nq
             Fmnu2(iq) = O1_NU(2)
             Fmnu0(iq) = O1_NU(0)
          enddo
          aniso = aniso +  pert%pa2_nu * sum(Fmnu2*wp)
       endif
       aniso = 0.6d0/ksq * aniso
       O1_PHI = O1_PSI - aniso

       !!velocities
       O1_V_C_PRIME = - O1_V_C + kbyaH * O1_PHI
       O1_NU_PRIME(1) = (O1_NU(0) + 4.d0*O1_PHI - 0.4d0 * O1_NU(2))*kbyaH
       do iq=1, pert%massivenu_iq_used
          O1_MASSIVENU_PRIME(1, iq) = (O1_MASSIVENU(0, iq)*qbye(iq) + 4.d0*O1_PHI/qbye(iq) - 0.4d0 * O1_MASSIVENU(2, iq)*qbye(iq)) * kbyaH
       enddo
       O1_V_B_PRIME = - O1_V_B + kbyaH * (O1_PHI + pert%cs2b * O1_DELTA_B) - pert%slip/(pert%R * aHtauc)
       O1_T_PRIME(1) = (O1_T(0) + 4.d0*O1_PHI - 0.4d0*pert%T%F(2))*kbyaH + 4.d0*pert%slip/aHtauc


       !!higher moments
       !!massless neutrinos
       do l = 2, pert%nu%lmax - 1
          O1_NU_PRIME(l) =  kbyaH * (cosmology%klms_by_2lm1(l, 0, 0) *   O1_NU( l-1 ) - cosmology%klms_by_2lp1(l+1, 0, 0) *  O1_NU( l+1 ) )
       enddo

       O1_NU_PRIME(pert%nu%lmax) = kbyaH *(pert%nu%lmax+0.5d0)/(pert%nu%lmax-0.5d0)*  O1_NU(pert%nu%lmax-1) &
            -  (pert%nu%lmax+1)* O1_NU(pert%nu%lmax)/(aHtau)

       if(cosmology%index_massivenu .ne. 0)then
          !!massive neutrinos
          do iq = 1, pert%massivenu_iq_used
             do l = 2, pert%massivenu(iq)%lmax - 1
                O1_MASSIVENU_PRIME(l, iq) = kbyaH * qbye(iq) * (cosmology%klms_by_2lm1(l, 0, 0) * O1_MASSIVENU(l-1, iq) - cosmology%klms_by_2lp1(l+1, 0, 0) * O1_MASSIVENU(l+1,iq))
             enddo
             O1_MASSIVENU_PRIME(pert%massivenu(iq)%lmax, iq) =  kbyaH * qbye(iq) * (pert%massivenu(iq)%lmax+0.5d0)/(pert%massivenu(iq)%lmax-0.5d0) *  O1_MASSIVENU(pert%nu%lmax-1, iq) &
                  -  (pert%nu%lmax+1)* O1_MASSIVENU(pert%nu%lmax, iq) / aHtau 
          enddo
       endif
       
       if(pert%tight_coupling)then

          pert%T2prime = (8.d0/9.d0)*(ktauc*O1_T_PRIME(1) + ktaucdot/pert%aH*O1_T(1)) !2.d0*pert%T%F(2)

          pert%E%F(2) = -coop_sqrt6/4.d0 * pert%T%F(2)
          pert%E2prime = -coop_sqrt6/4.d0 * pert%T2prime 

          aniso_prime =  pert%pa2_nu * O1_NU_PRIME(2) + pa2pr_nu*O1_NU(2) +  pert%pa2_g * pert%T2prime + pa2pr_g * pert%T%F(2)
       else
          !!T
          P = (O1_T(2) - coop_sqrt6 * O1_E(2))/10.d0
          O1_T_PRIME(2) =  kbyaH * (cosmology%klms_by_2lm1(2, 0, 0)*O1_T(1) - cosmology%klms_by_2lp1(3, 0, 0)*O1_T(3))  - (O1_T(2) - P)/aHtauc
          pert%T2prime =  O1_T_PRIME(2)
          do l = 3, pert%T%lmax -1
             O1_T_PRIME(l) = kbyaH * (cosmology%klms_by_2lm1(l, 0, 0)*O1_T(l-1) -cosmology%klms_by_2lp1(l+1, 0, 0)*O1_T(l+1))  - O1_T(l)/aHtauc
          enddo
          O1_T_PRIME(pert%T%lmax) =  kbyaH *((pert%T%lmax+0.5d0)/(pert%T%lmax-0.5d0))*  O1_T(pert%T%lmax-1) &
            -  ((pert%T%lmax+1)/aHtau + 1.d0/aHtauc) * O1_T(pert%T%lmax)
          !!E
          O1_E_PRIME(2) = kbyaH * ( - cosmology%klms_by_2lp1(3, 0, 0)*O1_E(3))  - (O1_E(2) + coop_sqrt6 * P)/aHtauc
          pert%E2prime = O1_E_PRIME(2)

          do l = 3, pert%E%lmax - 1
             O1_E_PRIME(l) = kbyaH * (cosmology%klms_by_2lm1(l, 0, 2)*O1_E(l-1) - cosmology%klms_by_2lp1(l+1, 0, 2)*O1_E(l+1))  - O1_E(l)/aHtauc
          enddo
          O1_E_PRIME(pert%E%lmax) =  kbyaH *( (pert%E%lmax+0.5d0)/(pert%E%lmax-0.5d0)) *  O1_E(pert%E%lmax-1) &
               -  ((pert%E%lmax+1)/aHtau + 1.d0/aHtauc) * O1_E(pert%E%lmax)
          aniso_prime = pert%pa2_g * O1_T_PRIME(2) + pa2pr_g * O1_T(2) + pert%pa2_nu * O1_NU_PRIME(2) + pa2pr_nu*O1_NU(2)
       endif
       if(cosmology%index_massivenu .ne. 0)then
          do iq = 1, pert%massivenu_iq_used
             Fmnu2_prime(iq) = O1_MASSIVENU_PRIME(2, iq)
          enddo
          do iq = pert%massivenu_iq_used + 1, coop_pert_default_nq
             Fmnu2_prime(iq) = O1_NU_PRIME(2)
          enddo
          aniso_prime = aniso_prime + pa2pr_nu * sum(Fmnu2*wp) + pert%pa2_nu*sum(Fmnu2_prime*wp + Fmnu2*wp_prime) 
       endif
       aniso_prime =  0.6d0/ksq * aniso_prime
       O1_PHI_PRIME = O1_PSI_PRIME - aniso_prime
       O1_PSIPR_PRIME = - O1_PHI_PRIME - (3.d0 + pert%daHdtau/aHsq)*O1_PSI_PRIME &
            - 2.d0*(pert%daHdtau/aHsq + 1.d0)*O1_PHI &
            - kbyaHsq/3.d0*(O1_PSI+aniso) &
            + (pert%rhoa2_b/aHsq * O1_DELTA_B * (pert%cs2b - 1.d0/3.d0) + pert%rhoa2_c/aHsq*O1_DELTA_C*(-1.d0/3.d0))/2.d0

       if(cosmology%index_massivenu .ne. 0)then
          pert%delta_mnu = sum(Fmnu0*wrho)
          pert%deltap_mnu = sum(Fmnu0*wp)
          O1_PSIPR_PRIME =  O1_PSIPR_PRIME + (pert%pa2_nu * pert%deltap_mnu - pert%rhoa2_nu/ 3.d0  * pert%delta_mnu )/2.d0
       else
          pert%delta_mnu = 0.d0
          pert%deltap_mnu = 0.d0
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
