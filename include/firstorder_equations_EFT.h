  subroutine coop_cosmology_firstorder_equations(n, lna, y, yp, cosmology, pert)
    COOP_INT n
    type(coop_cosmology_firstorder)::cosmology
    type(coop_pert_object)::pert
    COOP_REAL lna, y(0:n-1), yp(0:n-1)
    COOP_INT i, l, iq
    COOP_REAL a, aniso,  ktauc, ktaucdot, ktaucdd, aniso_prime, aHtauc, aHtau, aHsq, uterm, vterm, ma, doptdlna
    COOP_REAL :: pa2pr_g, pa2pr_nu
    COOP_REAL, dimension(coop_pert_default_nq)::Fmnu2_prime, Fmnu2, Fmnu0, qbye, wp, wp_prime, wrho_minus_wp
    COOP_REAL:: de_Q, de_dQdphi, de_phi,  de_V, de_Vp, de_Vpp, de_dphidlna,  asq, Hsq, de_m2byh2
    COOP_REAL,parameter::max_de_m2byH2 = 3.d2
    !!My PHI = Psi in Hu & White = Psi in Ma et al;
    !!My PSI = - Phi in Hu & White = Phi in Ma et al;
    !!My multipoles  = 4 * multipoles in Hu & White = (2l + 1) * multipoles in Ma et al
    !!My neutrino multipoles = (2l + 1) / (d ln f_0/d ln q) * Psi_l(q) in Ma et al
    !!time variable is log(a)
    yp(0) = 0
    pert%ksq = pert%k**2
    a = exp(lna)
    asq = a**2
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
    pert%M2 = cosmology%M2(a)
    pert%alpha_M = cosmology%alpha_M(a)
    pert%alpha_B = cosmology%alpha_B(a)
    pert%alpha_T = cosmology%alpha_T(a)
    pert%alpha_K = cosmology%alpha_K(a)
    pert%alpha_H = cosmology%alpha_H(a)

    pert%alpha_M_prime = cosmology%alpha_M_prime(a)
    pert%alpha_B_prime = cosmology%alpha_B_prime(a)
    pert%alpha_T_prime = cosmology%alpha_T_prime(a)
    pert%alpha_K_prime = cosmology%alpha_K_prime(a)
    pert%alpha_H_prime = cosmology%alpha_H_prime(a)

    
    if(cosmology%index_massivenu .ne. 0 )then
       pert%rhoa2_mnu = O0_MASSIVENU(cosmology)%rhoa2(a)
       pert%pa2_mnu = pert%rhoa2_mnu *  O0_MASSIVENU(cosmology)%wofa(a)
       ma = cosmology%mnu_by_Tnu * a
       if(ma  .lt. 1.d-4)then
          wrho_minus_wp = (cosmology%mnu_by_Tnu * a)**2/coop_pert_default_q**2
          qbye = 1.d0 - wrho_minus_wp/2.d0
          wrho_minus_wp = wrho_minus_wp * coop_pert_default_q_kernel *  pert%num_mnu_ratio
       else
          qbye =  coop_pert_default_q/sqrt(coop_pert_default_q**2 + (cosmology%mnu_by_Tnu * a)**2)
          wrho_minus_wp = coop_pert_default_q_kernel *  pert%num_mnu_ratio * (1.d0/qbye - qbye)
       endif
       wp = qbye * coop_pert_default_q_kernel * pert%num_mnu_ratio
       wp_prime = -(cosmology%mnu_by_Tnu*a)**2/(coop_pert_default_q**2 + (cosmology%mnu_by_Tnu * a)**2) * wp
    else
       pert%rhoa2_mnu = 0.d0
       pert%pa2_mnu = 0.d0
    endif

    pert%rhoa2_sum = pert%rhoa2_g + pert%rhoa2_b + pert%rhoa2_c + pert%rhoa2_nu +pert%rhoa2_mnu+pert%rhoa2_de 
    pert%pa2_sum = pert%pa2_g + pert%pa2_nu + pert%pa2_mnu + pert%pa2_de
    aHsq = (pert%rhoa2_sum + cosmology%Omega_k())/3.d0/pert%M2
    Hsq = aHsq/asq

    pert%aH = sqrt(aHsq)

    if(O0_DE(cosmology)%fDE_phi%initialized)then
       de_phi = O0_DE(cosmology)%DE_phi(a)
       pert%de_phidot = O0_DE(cosmology)%DE_phidot(a)
       de_Q = O0_DE(cosmology)%DE_Q(de_phi)
       de_dQdphi = O0_DE(cosmology)%DE_dlnQdphi(de_phi) * de_Q
       call  O0_DE(cosmology)%DE_get_VVpVpp(de_phi, de_V, de_Vp, de_Vpp)
       de_dphidlna = pert%de_phidot*a/pert%aH
       pert%de_Vp = de_Vp
    endif
    

    pert%daHdtau = -(pert%rhoa2_sum+3.d0*pert%pa2_sum)/6.d0/pert%M2
    pert%HdotbyHsq = pert%daHdtau/aHsq - 1.d0
    pert%HdotbyHsq_prime = cosmology%HddbyH3(a) - 2.d0*pert%HdotbyHsq**2    
    
    pa2pr_nu = O0_NU(cosmology)%dpa2da(a)*a
    pa2pr_g = O0_RADIATION(cosmology)%dpa2da(a)*a

    pert%R = 0.75d0 * pert%rhoa2_b/pert%rhoa2_g
    pert%tauc = cosmology%taucofa(a)
    pert%taucdot = cosmology%dot_tauc(a)

    ktauc = pert%k * pert%tauc
    ktaucdot = pert%k * pert%taucdot
    pert%kbyaH  = pert%k/pert%aH
    pert%kbyaHsq = pert%kbyaH**2
    aHtauc = pert%aH * pert%tauc
    doptdlna = 1.d0/aHtauc
    pert%tau = cosmology%tauofa(a)
    aHtau = pert%aH*pert%tau

    ktaucdd = pert%k*(cosmology%dot_tauc(a*1.01)-cosmology%dot_tauc(a/1.01))/(0.02/pert%aH)					      
    select case(pert%m)
    case(0)
       O1_PSI_PRIME = O1_PSIPR
       O1_DELTA_C_PRIME = - O1_V_C * pert%kbyaH + 3.d0 * O1_PSI_PRIME
       O1_DELTA_B_PRIME = - O1_V_B * pert%kbyaH + 3.d0 * O1_PSI_PRIME
       if(pert%has_rad_pert) &
            O1_T_PRIME(0) = - O1_T(1) * pert%kbyaH/3.d0 + 4.d0 * O1_PSI_PRIME 
       O1_NU_PRIME(0) = - O1_NU(1) * pert%kbyaH/3.d0 + 4.d0 * O1_PSI_PRIME

       if(pert%tight_coupling)then
          pert%T%F(2) = (8.d0/9.d0)*ktauc * O1_T(1)
          Uterm = - pert%aH * O1_V_B + (pert%cs2b * O1_DELTA_B - O1_T(0)/4.d0 +  pert%T%F(2)/10.d0)*pert%k
          vterm = pert%k * (pert%R+1.d0)/pert%R + ktaucdot 
          pert%slip = ktauc &
               * (Uterm - ktauc*(Uterm*(-pert%aH*pert%k/pert%R + pert%daHdtau/pert%aH*ktaucdot)/vterm)/vterm)/vterm   !!v_b - v_g accurate to (k tau_c)^2 
       else
          if(pert%has_rad_pert)then          
             pert%T%F(2) = O1_T(2)
             pert%slip = O1_V_B - O1_T(1)/4.d0
          else
             pert%T%F(2) = 0.d0
             pert%slip = 0.d0
          endif
       endif
       aniso = (pert%pa2_g * pert%T%F(2) + pert%pa2_nu * O1_NU(2))

       if(cosmology%index_massivenu .ne. 0 .and. .not. pert%massivenu_cold)then
          do iq = 1, pert%massivenu_iq_used
             Fmnu2(iq) = O1_MASSIVENU(2, iq)
             Fmnu0(iq) = O1_MASSIVENU(0, iq)
             O1_MASSIVENU_PRIME(0, iq) = - O1_MASSIVENU(1, iq)*pert%kbyaH/3.d0 * qbye(iq) + 4.d0 * O1_PSI_PRIME 
          enddo
          do iq = pert%massivenu_iq_used + 1, coop_pert_default_nq
             Fmnu2(iq) = O1_NU(2)
             Fmnu0(iq) = O1_NU(0)
          enddo
          aniso = aniso +  pert%pa2_nu * sum(Fmnu2*wp)
       endif
       aniso = 0.6d0/pert%ksq * aniso
       O1_PHI = (O1_PSI*(1.d0+pert%alpha_T) - aniso/pert%M2 - (pert%alpha_M - pert%alpha_T) * O1_DE_HPI + pert%alpha_H * (O1_DE_HPI_PRIME - pert%HdotbyHsq * O1_DE_HPI))/(1.d0+pert%alpha_H)

       !!velocities
       O1_V_C_PRIME = - O1_V_C + pert%kbyaH * O1_PHI
       O1_NU_PRIME(1) = (O1_NU(0) + 4.d0*O1_PHI - 0.4d0 * O1_NU(2))*pert%kbyaH
       do iq=1, pert%massivenu_iq_used
          O1_MASSIVENU_PRIME(1, iq) = (O1_MASSIVENU(0, iq)*qbye(iq) + 4.d0*O1_PHI/qbye(iq) - 0.4d0 * O1_MASSIVENU(2, iq)*qbye(iq)) * pert%kbyaH
       enddo
       O1_V_B_PRIME = - O1_V_B + pert%kbyaH * (O1_PHI + pert%cs2b * O1_DELTA_B) - pert%slip/(pert%R * aHtauc)
       if(pert%has_rad_pert) &       
            O1_T_PRIME(1) = (O1_T(0) + 4.d0*O1_PHI - 0.4d0*pert%T%F(2))*pert%kbyaH + 4.d0*pert%slip/aHtauc


       !!higher moments
       !!massless neutrinos
       do l = 2, pert%nu%lmax - 1
          O1_NU_PRIME(l) =  pert%kbyaH * (cosmology%klms_by_2lm1(l, 0, 0) *   O1_NU( l-1 ) - cosmology%klms_by_2lp1(l+1, 0, 0) *  O1_NU( l+1 ) )
       enddo

       O1_NU_PRIME(pert%nu%lmax) = pert%kbyaH * (pert%nu%lmax +0.5d0)/(pert%nu%lmax-0.5d0)*  O1_NU(pert%nu%lmax-1) -  (pert%nu%lmax+1)* O1_NU(pert%nu%lmax)/(aHtau)


       if(cosmology%index_massivenu .ne. 0)then
          !!massive neutrinos
          if(pert%massivenu_cold)then
             O1_MASSIVENU_PRIME(0, 1) = - O1_MASSIVENU(1, 1)* pert%kbyaH+ 3.d0 * O1_PSI_PRIME
             O1_MASSIVENU_PRIME(1, 1) = - O1_MASSIVENU(1, 1) + pert%kbyaH * O1_PHI
          else
             do iq = 1, pert%massivenu_iq_used
                do l = 2, pert%massivenu(iq)%lmax - 1
                   O1_MASSIVENU_PRIME(l, iq) = pert%kbyaH * qbye(iq) * (cosmology%klms_by_2lm1(l, 0, 0) * O1_MASSIVENU(l-1, iq) - cosmology%klms_by_2lp1(l+1, 0, 0) * O1_MASSIVENU(l+1,iq))
                enddo
                O1_MASSIVENU_PRIME(pert%massivenu(iq)%lmax, iq) =  pert%kbyaH * qbye(iq) * (pert%massivenu(iq)%lmax+0.5d0)/(pert%massivenu(iq)%lmax-0.5d0) *  O1_MASSIVENU(pert%nu%lmax-1, iq) &
                     -  (pert%nu%lmax+1)* O1_MASSIVENU(pert%nu%lmax, iq) / aHtau 
             enddo
          endif
       endif
       
       if(pert%tight_coupling)then

          pert%T2prime = (8.d0/9.d0)*(ktauc*O1_T_PRIME(1) + ktaucdot/pert%aH*O1_T(1)) !2.d0*pert%T%F(2)
          pert%E%F(2) = -coop_sqrt6/4.d0 * pert%T%F(2)
          pert%E2prime = -coop_sqrt6/4.d0 * pert%T2prime 
          pert%capP = (pert%T%F(2) - coop_sqrt6 * pert%E%F(2))/10.d0

       else
          !!T
          if(pert%has_rad_pert)then
             pert%capP = (O1_T(2) - coop_sqrt6 * O1_E(2))/10.d0
             O1_T_PRIME(2) =  pert%kbyaH * (cosmology%klms_by_2lm1(2, 0, 0)*O1_T(1) - cosmology%klms_by_2lp1(3, 0, 0)*O1_T(3))  - (O1_T(2) - pert%capP)/aHtauc
             pert%T2prime =  O1_T_PRIME(2)
             do l = 3, pert%T%lmax -1
                O1_T_PRIME(l) = pert%kbyaH * (cosmology%klms_by_2lm1(l, 0, 0)*O1_T(l-1) -cosmology%klms_by_2lp1(l+1, 0, 0)*O1_T(l+1))  - O1_T(l)/aHtauc
             enddo
             O1_T_PRIME(pert%T%lmax) =  pert%kbyaH *((pert%T%lmax+0.5d0)/(pert%T%lmax-0.5d0))*  O1_T(pert%T%lmax-1) &
                  -  ((pert%T%lmax+1)/aHtau + doptdlna) * O1_T(pert%T%lmax)
          !!E
             O1_E_PRIME(2) = pert%kbyaH * ( - cosmology%klms_by_2lp1(3, 0, 2)*O1_E(3))  - (O1_E(2) + coop_sqrt6 * pert%capP)/aHtauc
             pert%E2prime = O1_E_PRIME(2)

             do l = 3, pert%E%lmax - 1
                O1_E_PRIME(l) = pert%kbyaH * (cosmology%klms_by_2lm1(l, 0, 2)*O1_E(l-1) - cosmology%klms_by_2lp1(l+1, 0, 2)*O1_E(l+1)) - O1_E(l)/aHtauc
             enddo
             O1_E_PRIME(pert%E%lmax) =  pert%kbyaH *( (pert%E%lmax-2.d0/pert%E%lmax+0.5d0)/(pert%E%lmax-2.d0/pert%E%lmax-0.5d0)) *  O1_E(pert%E%lmax-1) &
                  -  ((pert%E%lmax-2.d0/pert%E%lmax+1)/aHtau + doptdlna) * O1_E(pert%E%lmax)

          else
             pert%capP = 0.d0
             pert%T2prime = 0.d0
             pert%E2prime = 0.d0
          endif
       endif
       aniso_prime =  pert%pa2_nu * O1_NU_PRIME(2) + pa2pr_nu*O1_NU(2) +  pert%pa2_g * pert%T2prime + pa2pr_g * pert%T%F(2)       
       if(cosmology%index_massivenu .ne. 0 .and. .not. pert%massivenu_cold)then
          do iq = 1, pert%massivenu_iq_used
             Fmnu2_prime(iq) = O1_MASSIVENU_PRIME(2, iq)
          enddo
          do iq = pert%massivenu_iq_used + 1, coop_pert_default_nq
             Fmnu2_prime(iq) = O1_NU_PRIME(2)
          enddo
          aniso_prime = aniso_prime + pa2pr_nu * sum(Fmnu2*wp) + pert%pa2_nu*sum(Fmnu2_prime*wp + Fmnu2*wp_prime) 
       endif
       aniso_prime =  0.6d0/pert%ksq * aniso_prime
       
!!$       O1_PHI_PRIME = O1_PSI_PRIME - aniso_prime
       O1_PHI_PRIME = O1_PHI*(-pert%alpha_H_prime)/(1.d0+pert%alpha_H) &
            + ( &
            O1_PSIPR*(1.d0+pert%alpha_T) + O1_PSI*(1.d0+pert%alpha_T_prime) &
            - aniso_prime/pert%M2  - aniso/pert%M2*pert%alpha_M &
            - (pert%alpha_M_prime - pert%alpha_T_prime) * O1_DE_HPI - (pert%alpha_M - pert%alpha_T) * O1_DE_HPIPR &
            + pert%alpha_H * ( - pert%HdotbyHsq * O1_DE_HPIPR - pert%HdotbyHsq_prime*O1_DE_HPI) + pert%alpha_H_prime * (O1_DE_HPI_PRIME - pert%HdotbyHsq * O1_DE_HPI) &
            )/(1.d0+pert%alpha_H)  !!this is actually Phi_Prime - alpha_H/(1+alpha_H) PsiPr_Prime
       
       if(pert%want_source)then
          pert%ekappa = cosmology%ekappaofa(pert%a)
          pert%vis = pert%ekappa/pert%tauc
          pert%visdot = cosmology%vis%derivative(pert%a) * pert%a * pert%aH
          pert%Pdot = (pert%T2prime  - coop_sqrt6 *  pert%E2prime)/10.d0 * pert%aH
          pert%kchi = 1.d0 - pert%tau/cosmology%tau0
          pert%kchi = pert%k * cosmology%tau0 * (pert%kchi + exp(-1.d3*pert%kchi))
       endif
       O1_PSIPR_PRIME = - (1.d0+pert%alpha_B)*O1_PHI_PRIME &
            - (3.d0 + pert%alpha_M  + pert%daHdtau/aHsq)*O1_PSI_PRIME &
            - 2.d0*(pert%daHdtau/aHsq + 1.d0)*O1_PHI &
            - pert%kbyaHsq/3.d0*(O1_PSI+aniso/pert%M2) &
            + ( &
            pert%rhoa2_b/aHsq * O1_DELTA_B * (pert%cs2b - 1.d0/3.d0) &
            + pert%rhoa2_c/aHsq*O1_DELTA_C*(-1.d0/3.d0) &
            )/2.d0/pert%M2

       if(cosmology%index_massivenu .ne. 0)then
          if(pert%massivenu_cold)then
             O1_PSIPR_PRIME = O1_PSIPR_PRIME - O1_MASSIVENU(0, 1)*pert%rhoa2_mnu/6.0/aHsq
          else
             pert%deltatr_mnu = sum(Fmnu0*wrho_minus_wp)
	     pert%deltap_mnu = sum(Fmnu0*wp)
             O1_PSIPR_PRIME =  O1_PSIPR_PRIME - (pert%pa2_nu/aHsq * pert%deltatr_mnu)/2.d0
          endif
       else
          pert%deltatr_mnu = 0.d0
          pert%deltap_mnu = 0.d0
       endif

       select case(pert%de%genre)
       case(COOP_PERT_NONE)
          !!do nothing
       case(COOP_PERT_SCALAR_FIELD)  !!coupled DE
          O1_DELTA_C_PRIME =   O1_DELTA_C_PRIME  + de_Q * O1_DELTA_PHIPR + de_dQdphi*O1_DELTA_PHI * de_dphidlna
          O1_V_C_PRIME =   O1_V_C_PRIME - de_Q * de_dphidlna * O1_V_C + pert%kbyaH * de_Q * O1_DELTA_PHI
          O1_DELTA_PHI_PRIME = O1_DELTA_PHIPR
          de_m2byH2 = pert%kbyaHsq + de_Vpp/Hsq
          if(de_m2byH2 .lt. max_de_m2byH2)then
             O1_DELTA_PHIPR_PRIME =  - (2.d0 + pert%daHdtau/aHsq)* O1_DELTA_PHIPR -de_m2byH2 *O1_DELTA_PHI + (3.d0*O1_PSIPR + O1_PHI_PRIME)*de_dphidlna - 2.d0*de_Vp/Hsq*O1_PHI - pert%rhoa2_c/aHsq* (de_Q*(O1_DELTA_C + 2.d0*O1_PHI) + de_dQdphi*O1_DELTA_PHI)
          else
             O1_DELTA_PHIPR_PRIME =  - (2.d0 + pert%daHdtau/aHsq)* O1_DELTA_PHIPR - max_de_m2byH2 *O1_DELTA_PHI + ((3.d0*O1_PSIPR + O1_PHI_PRIME)*de_dphidlna - 2.d0*de_Vp/Hsq*O1_PHI - pert%rhoa2_c/aHsq* (de_Q*(O1_DELTA_C + 2.d0*O1_PHI) + de_dQdphi*O1_DELTA_PHI))*max_de_m2byH2/de_m2byH2             
          endif
          pert%de_delta_rho = pert%de_phidot * pert%O1_DELTA_PHIPR * pert%aH / pert%a &
               - (pert%de_phidot)**2*O1_PHI + pert%de_Vp * O1_DELTA_PHI
          pert%de_delta_p = pert%de_phidot * pert%O1_DELTA_PHIPR * pert%aH / pert%a &
               - (pert%de_phidot)**2*O1_PHI - pert%de_Vp * O1_DELTA_PHI
          O1_PSIPR_PRIME = O1_PSIPR_PRIME  + ( &
               (pert%de_delta_rho*(-1.d0/3.d0)+pert%de_delta_p)/Hsq &
            )/2.d0
       case(COOP_PERT_EFT)
          O1_DE_HPI_PRIME = O1_DE_HPIPR
          O1_DE_HPIPR_PRIME = 0.d0
       case default
          call coop_tbw("DE perturbations not written")
       end select
    case(1)
       call coop_tbw("vector equations not written")
    case(2)
       O1_TEN_H_PRIME = O1_TEN_HPR
       if(pert%tight_coupling)then
          pert%T%F(2) =  (-16.d0/3.d0 )*O1_TEN_HPR*aHtauc
!(-16.d0/3.d0 + (8.d0*19.d0/63.d0)*(pert%kbyaH**2*aHtauc + (231.D0/19.D0))*aHtauc  )*O1_TEN_HPR*aHtauc
          pert%E%F(2) =  (-coop_sqrt6/4.d0)*pert%T%F(2) !- (110.d0/63.d0*coop_sqrt6)*(pert%kbyaH**2*aHtauc +(63.d0/11.d0))*O1_TEN_HPR*aHtauc**2
       else
          if(pert%has_rad_pert)then
             pert%T%F(2) =O1_T(2)
             pert%E%F(2) = O1_E(2)
          endif
       endif
       if(pert%has_rad_pert)then       
          pert%capP = (pert%T%F(2) - coop_sqrt6 * pert%E%F(2))/10.d0
          aniso = pert%pa2_g * pert%T%F(2) + pert%pa2_nu * O1_NU(2)          
       else
          pert%capP = 0.d0
          aniso =  pert%pa2_nu * O1_NU(2)          
       endif
       if(cosmology%index_massivenu .ne. 0 .and. .not. pert%massivenu_cold)then
          do iq = 1, pert%massivenu_iq_used
             Fmnu2(iq) = O1_MASSIVENU(2, iq)
          enddo
          do iq = pert%massivenu_iq_used + 1, coop_pert_default_nq
             Fmnu2(iq) = O1_NU(2)
          enddo
          aniso = aniso +  pert%pa2_nu * sum(Fmnu2*wp)
       endif
       O1_TEN_HPR_PRIME = -(2.d0+pert%daHdtau/aHsq)*O1_TEN_HPR &
            - pert%kbyaH**2 * O1_TEN_H &
            + 0.4d0/aHsq * aniso


       O1_NU_PRIME(2) =  pert%kbyaH * ( - cosmology%klms_by_2lp1(3, 2, 0) *  O1_NU( 3 ) ) - 4.d0*O1_TEN_HPR
       do l = 3, pert%nu%lmax - 1
          O1_NU_PRIME(l) =  pert%kbyaH * (cosmology%klms_by_2lm1(l, 2, 0) *   O1_NU( l-1 ) - cosmology%klms_by_2lp1(l+1, 2, 0) *  O1_NU( l+1 ) )
       enddo

       O1_NU_PRIME(pert%nu%lmax) = pert%kbyaH * (pert%nu%lmax-2.d0/pert%nu%lmax+0.5d0)/(pert%nu%lmax-2.d0/pert%nu%lmax-0.5d0)*  O1_NU(pert%nu%lmax-1) &
            -  (pert%nu%lmax-2.d0/pert%nu%lmax+1)* O1_NU(pert%nu%lmax)/(aHtau)

       if(cosmology%index_massivenu .ne. 0 .and. .not. pert%massivenu_cold)then
          !!massive neutrinos
          do iq = 1, pert%massivenu_iq_used
             O1_MASSIVENU_PRIME(2, iq) = pert%kbyaH * qbye(iq) * (- cosmology%klms_by_2lp1(3, 2, 0) * O1_MASSIVENU(3, iq))-4.d0*O1_TEN_HPR
             do l = 3, pert%massivenu(iq)%lmax - 1
                O1_MASSIVENU_PRIME(l, iq) = pert%kbyaH * qbye(iq) * (cosmology%klms_by_2lm1(l, 2, 0) * O1_MASSIVENU(l-1, iq) - cosmology%klms_by_2lp1(l+1, 2, 0) * O1_MASSIVENU(l+1,iq))
             enddo
             O1_MASSIVENU_PRIME(pert%massivenu(iq)%lmax, iq) =  pert%kbyaH * qbye(iq) * (pert%massivenu(iq)%lmax-2.d0/pert%massivenu(iq)%lmax+0.5d0)/(pert%massivenu(iq)%lmax-2.d0/pert%massivenu(iq)%lmax-0.5d0) *  O1_MASSIVENU(pert%nu%lmax-1, iq) &
                  -  (pert%nu%lmax-2.d0/pert%massivenu(iq)%lmax+1)* O1_MASSIVENU(pert%nu%lmax, iq) / aHtau 
          enddo
       endif

       if(pert%tight_coupling)then
          pert%T2prime = -16.d0/3.d0*(O1_TEN_HPR_PRIME*aHtauc - O1_TEN_HPR * (pert%taucdot  + pert%tauc*pert%daHdtau/pert%aH)/aHtauc**2)
          pert%E2prime = (-coop_sqrt6/4.d0)*pert%T2prime
       else
          if(pert%has_rad_pert)then
             O1_T_PRIME(2) =  pert%kbyah * ( -  cosmology%klms_by_2lp1(3, 2, 0) * O1_T(3))  -  ( O1_T(2) - pert%capP)/aHtauc &
                  - O1_TEN_HPR*4.d0
             pert%T2prime = O1_T_PRIME(2)
             O1_E_PRIME(2) = pert%kbyah * ( - cosmology%fourbyllp1(2)*O1_B(2)- cosmology%klms_by_2lp1(3, 2, 2) * O1_E(3)) - (O1_E(2) + coop_sqrt6 * pert%capP)/aHtauc
             pert%E2prime = O1_E_PRIME(2)
             O1_B_PRIME(2) = pert%kbyah * ( cosmology%fourbyllp1(2)*O1_E(2)- cosmology%klms_by_2lp1(3, 2, 2) * O1_B(3)) - O1_B(2)/aHtauc
             do l = 3, pert%T%lmax -1 
                O1_T_PRIME(l) = pert%kbyah * (cosmology%klms_by_2lm1(l, 2, 0) * O1_T(l-1) -  cosmology%klms_by_2lp1(l+1, 2, 0) * O1_T(l+1)) -  O1_T(l)/aHtauc
             enddo
             O1_T_PRIME(pert%T%lmax) = pert%kbyah  *((pert%T%lmax-2.d0/pert%T%lmax+0.5d0)/(pert%T%lmax-2.d0/pert%T%lmax-0.5d0)) * O1_T(pert%T%lmax-1) &
                  - (doptdlna + (pert%T%lmax-2.d0/pert%T%lmax+1)/aHtau)* O1_T(pert%T%lmax)
             if(pert%E%lmax .ne. pert%B%lmax) stop "lmax(E) must be the same as lmax(B)"
             do l=3, pert%E%lmax -1
                O1_E_PRIME(l) = pert%kbyah * (cosmology%klms_by_2lm1(l, 2, 2) * O1_E(l-1) &
                     - cosmology%fourbyllp1(l)*O1_B(l) &
                     - cosmology%klms_by_2lp1(l+1, 2, 2) * O1_E(l+1)) - O1_E(l)/aHtauc
                O1_B_PRIME(l) = pert%kbyah * (cosmology%klms_by_2lm1(l, 2, 2) * O1_B(l-1) + cosmology%fourbyllp1(l)*O1_E(l)- cosmology%klms_by_2lp1(l+1, 2, 2) * O1_B(l+1)) -  O1_B(l)/aHtauc
             enddo
             O1_E_PRIME(pert%E%lmax) = pert%kbyah  *((pert%E%lmax-4.d0/pert%E%lmax+0.5d0)/(pert%E%lmax-4.d0/pert%E%lmax-0.5d0)) * O1_E(pert%E%lmax-1) &
                  - cosmology%fourbyllp1(pert%E%lmax)*O1_B(pert%E%lmax) &
                  - (doptdlna + (pert%E%lmax-4.d0/pert%E%lmax+1)/aHtau)* O1_E(pert%E%lmax)
             O1_B_PRIME(pert%E%lmax) = pert%kbyah  *((pert%E%lmax-4.d0/pert%E%lmax+0.5d0)/(pert%E%lmax-4.d0/pert%E%lmax-0.5d0)) * O1_B(pert%E%lmax-1)  + cosmology%fourbyllp1(pert%E%lmax)*O1_E(pert%E%lmax) &
                  - (doptdlna + (pert%E%lmax-4.d0/pert%E%lmax+1)/aHtau)* O1_B(pert%E%lmax)
          else
             pert%T2prime = 0.d0
             pert%E2prime = 0.d0
          endif
       endif
    case default
       call coop_return_error("firstorder_equations", "Unknown m = "//trim(coop_num2str(pert%m)), "stop")
    end select
  end subroutine coop_cosmology_firstorder_equations