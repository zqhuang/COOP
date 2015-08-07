  subroutine coop_cosmology_firstorder_equations(n, lna, y, yp, cosmology, pert)
    COOP_INT n
    type(coop_cosmology_firstorder)::cosmology
    type(coop_pert_object)::pert
    COOP_REAL, parameter::tol = 1.d-5, conv_slope = 50.d0
    COOP_REAL lna, y(0:n-1), yp(0:n-1)
    COOP_INT i, l, iq
    COOP_REAL a, aniso,  ktauc, ktaucdot, ktaucdd, aniso_prime, aHtauc, aHtau, aHsq, uterm, vterm, ma, doptdlna, M2a2H2, anisobyM2
    COOP_REAL :: pa2pr_g, pa2pr_nu
    COOP_REAL, dimension(coop_pert_default_nq)::Fmnu2_prime, Fmnu2, Fmnu0, qbye, wp, wp_prime, wrho_minus_wp
    COOP_REAL::  asq, Hsq
    COOP_INT::category
    COOP_REAL:: Pc, P1, P0, Pcpr, P1pr, P0pr, phi_c, phipr_c, pidot_c, Mat(4, 2), cons(2)
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

    pert%rhoa2_mnu = 0.d0
    pert%pa2_mnu = 0.d0

    pert%rhoa2_matter = pert%rhoa2_g + pert%rhoa2_b + pert%rhoa2_c + pert%rhoa2_nu +pert%rhoa2_mnu
    pert%pa2_matter = pert%pa2_g + pert%pa2_nu + pert%pa2_mnu

    
    pert%rhoa2_sum = pert%rhoa2_g + pert%rhoa2_b + pert%rhoa2_c + pert%rhoa2_nu +pert%rhoa2_mnu+pert%rhoa2_de 
    pert%pa2_sum = pert%pa2_g + pert%pa2_nu + pert%pa2_mnu + pert%pa2_de
    aHsq = (pert%rhoa2_sum + cosmology%Omega_k())/3.d0/pert%M2  !!a^2H^2
    M2a2H2 = pert%M2 * aHsq
    Hsq = aHsq/asq

    pert%aH = sqrt(aHsq)

    pert%daHdtau = -(pert%rhoa2_sum+3.d0*pert%pa2_sum)/6.d0/pert%M2
    pert%HdotbyHsq = pert%daHdtau/aHsq - 1.d0
    pert%HddbyH3 = cosmology%HddbyH3(a)
    pert%HdotbyHsq_prime = pert%HddbyH3 - 2.d0*pert%HdotbyHsq**2    
    
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
       aniso = 0.6d0/pert%ksq * aniso
       anisobyM2 = aniso/pert%M2


       !!Phi = Pc + P0 Hpi + P1 Hpi'
       Pc = ((1.d0+pert%alpha_T)*O1_PSI - anisobyM2)/(1.d0+pert%alpha_H)
       P1 = pert%alpha_H/(1.d0+pert%alpha_H)
       P0 = (pert%alpha_T - pert%alpha_M - pert%HdotbyHsq*pert%alpha_H)/(1.d0+pert%alpha_H)
       
       !!Eq. (111) 
       O1_PHI = Pc + P0 *  O1_DE_HPI + P1 * O1_DE_HPIPR

       !!velocities
       O1_V_C_PRIME = - O1_V_C + pert%kbyaH * O1_PHI
       O1_NU_PRIME(1) = (O1_NU(0) + 4.d0*O1_PHI - 0.4d0 * O1_NU(2))*pert%kbyaH
       O1_V_B_PRIME = - O1_V_B + pert%kbyaH * (O1_PHI + pert%cs2b * O1_DELTA_B) - pert%slip/(pert%R * aHtauc)
       if(pert%has_rad_pert) &       
            O1_T_PRIME(1) = (O1_T(0) + 4.d0*O1_PHI - 0.4d0*pert%T%F(2))*pert%kbyaH + 4.d0*pert%slip/aHtauc


       !!higher moments
       !!massless neutrinos
       do l = 2, pert%nu%lmax - 1
          O1_NU_PRIME(l) =  pert%kbyaH * (cosmology%klms_by_2lm1(l, 0, 0) *   O1_NU( l-1 ) - cosmology%klms_by_2lp1(l+1, 0, 0) *  O1_NU( l+1 ) )
       enddo

       O1_NU_PRIME(pert%nu%lmax) = pert%kbyaH * (pert%nu%lmax +0.5d0)/(pert%nu%lmax-0.5d0)*  O1_NU(pert%nu%lmax-1) -  (pert%nu%lmax+1)* O1_NU(pert%nu%lmax)/(aHtau)


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
       aniso_prime = pert%pa2_nu * O1_NU_PRIME(2) + pa2pr_nu*O1_NU(2) +  pert%pa2_g * pert%T2prime + pa2pr_g * pert%T%F(2)
       aniso_prime =  0.6d0/pert%ksq * aniso_prime

       !!now start doing the EFT part
       Pcpr = (- pert%alpha_H_prime*Pc + pert%alpha_T_prime * O1_PSI + (1.d0+pert%alpha_T)*O1_PSIPR + (aniso * pert%alpha_M - aniso_prime)/pert%M2)/(1.d0+pert%alpha_H)
       P1pr = pert%alpha_H_prime/(1.d0+pert%alpha_H)**2
       P0pr = (-pert%alpha_H_prime*P0 + pert%alpha_T_prime - pert%alpha_M_prime - pert%HdotbyHsq * pert%alpha_H_prime - pert%HdotbyHsq_prime*pert%alpha_H)/(1.d0+pert%alpha_H)

       pert%u = - (pert%rhoa2_de + pert%pa2_de)/(2.d0*M2a2h2)       

       Mat(1, 1) = 1.d0
       
       phi_c = 1.d0+2.d0*pert%HdotbyHsq - pert%u - pert%alpha_K/6.d0 + (2.d0+pert%HdotbyHsq)*pert%alpha_B + pert%alpha_B_prime + (3.d0+pert%alpha_M)*(1.d0+pert%alpha_B)
       phipr_c = 1.d0 + pert%alpha_B
       pidot_c = pert%u - pert%alpha_B_prime - (3.d0+pert%alpha_M + pert%HdotbyHsq)*pert%alpha_B
       Mat(2, 1) = P1 *  phipr_c - pert%alpha_B       
       Mat(3,1) = P1*phi_c + (P1pr + P0)*phipr_c + pert%alpha_B * pert%HdotbyHsq &
            + pidot_c + (pert%alpha_K/6.d0 - pert%alpha_B)
       Mat(4,1) = phi_c * P0 + phipr_c * P0pr &
            + pert%alpha_B * pert%HdotbyHsq_prime &
            - pert%HdotbyHsq * pidot_c &
            + (3.d0+pert%alpha_M + 2.d0*pert%hdotbyHsq)*pert%hdotbyHsq - 2.d0*pert%pa2_matter/M2a2H2 + pert%HdotbyHsq_prime  &
            - (pert%alpha_K/6.d0 - pert%alpha_B) * pert%HdotbyHsq  &
            + pert%u + pert%HdotbyHsq * pert%alpha_B + pert%kbyaHsq*(pert%alpha_H - pert%alpha_B)

       Cons(1) =  - phi_c * Pc - phipr_c * Pcpr + 0.5d0 * (O1_DELTA_B*pert%rhoa2_b/M2a2H2*(pert%cs2b - 1.d0/3.d0) + O1_DELTA_C * pert%rhoa2_c/M2a2H2*(-1.d0/3.d0)- 2.d0/3.d0 * pert%kbyaHsq * anisobyM2) &
               - (4.d0 + pert%HdotbyHsq + pert%alpha_M + pert%alpha_B)*O1_PSIPR &
               -((1.d0 + pert%alpha_H)/3.d0*pert%kbyaHsq)*O1_PSI

       phi_c  = 6.d0*pert%u + (6.d0*pert%alpha_B - pert%alpha_K) * (3.d0 + pert%alpha_M) + 2.d0*(pert%alpha_B * 9.d0 - pert%alpha_K) * pert%HdotbyHsq  + 6.d0*pert%alpha_B_prime - pert%alpha_K_prime + 2.d0*pert%kbyaHsq*(pert%alpha_H - pert%alpha_B)
       phipr_c = 6.d0*pert%alpha_B- pert%alpha_K
       Mat(1, 2) = 6.d0*pert%alpha_B
       Cons(2) = phi_c * Pc   + phipr_c * Pcpr &
            - (6.d0*(pert%alpha_B*(2.d0*pert%HdotbyHsq+3.d0+pert%alpha_M) + pert%u + pert%alpha_B_prime)+2.d0*pert%kbyaHsq*pert%alpha_H)*O1_PSIPR &
            - (2.d0*pert%kbyaHsq * (pert%alpha_M + pert%alpha_H + pert%alpha_M * pert%alpha_H - pert%alpha_T - pert%alpha_H_prime))*O1_PSI
       Mat(2, 2) = pert%alpha_K + phipr_c * P1
       Mat(3, 2) = - pert%alpha_K * pert%HdotbyHsq &
            + ((3.d0+pert%alpha_M + 2.d0*pert%HdotbyHsq)*pert%alpha_K + pert%alpha_K_prime) &
            + phipr_c * (P0 + P1pr) + phi_c * P1
       Mat(4, 2) = -pert%alpha_K * pert%HdotbyHsq_prime - pert%HdotbyHsq * ((3.d0 + pert%alpha_M + 2.d0*pert%HdotbyHsq)*pert%alpha_K + pert%alpha_K_prime ) &
            + 6.d0*(pert%u * pert%HdotbyHsq + pert%alpha_B * pert%HdotbyHsq * (3.d0 + pert%alpha_M  + pert%HdotbyHsq) + pert%HdotbyHsq * pert%alpha_B_prime + pert%alpha_B*(pert%HdotbyHsq_prime + 2.d0 * pert%HdotbyHsq ** 2) ) &
            - 2.d0 * pert%kbyaHsq * ( pert%u + 1.d0 + pert%alpha_B*(1.d0+pert%alpha_M) +pert%alpha_T -(1.d0+pert%alpha_H)*(1.d0+pert%alpha_M) + pert%HdotbyHsq*(pert%alpha_B - pert%alpha_H) + pert%alpha_B_prime - pert%alpha_H_prime) &
            + phipr_c * P0pr &
            + phi_c * P0

       
       Mat(:, 2) = Mat(:, 2) - Mat(:, 1)* Mat(1, 2) !!this makes Mat(1, 2) = 0
       Cons(2) = Cons(2) - Cons(1) * Mat(1, 2)
       if(abs(Mat(2, 2)) .lt. tol)then
          if( abs(Mat(3, 2)) .lt. tol)then
             if(abs(Mat(4,2)).lt. tol)then  !!no dof, LCDM case                
                O1_DE_HPIPR_PRIME = 0.d0
                O1_DE_HPI_PRIME = 0.d0
             else
                O1_DE_HPI_PRIME = - (O1_DE_HPI - Cons(2)/Mat(4,2))*conv_slope
                O1_DE_HPIPR_PRIME = -(O1_DE_HPIPR -  &
                     (Cons(1)  - Mat(3,1)*O1_DE_HPIPR - Mat(4,1)*O1_DE_HPI + PcPr  + (2.d0/3.d0*pert%kbyaHsq*(1.d0+pert%HdotbyHsq) + pert%HdotbyHsq_prime)/(pert%kbyaHsq/3.d0- pert%HdotbyHsq)*(O1_PSIPR + O1_PHI))/(pert%kbyaHsq/3.d0- pert%HdotbyHsq) &
                     )*conv_slope
             endif
          else
             O1_DE_HPI_PRIME = O1_DE_HPIPR
             O1_DE_HPIPR_PRIME = - (O1_DE_HPIPR - (Cons(2) - O1_DE_HPI*Mat(4,2))/Mat(3, 2))*conv_slope
          endif
       else
          O1_DE_HPI_PRIME = O1_DE_HPIPR
          O1_DE_HPIPR_PRIME = (Cons(2) - O1_DE_HPIPR*Mat(3, 2) - O1_DE_HPI*Mat(4, 2))/Mat(2, 2)
       endif
       O1_PSIPR_PRIME = Cons(1) - Mat(2,1)*O1_DE_HPIPR_PRIME - Mat(3,1)*O1_DE_HPIPR - Mat(4,1)*O1_DE_HPI
       if(pert%want_source)then
          pert%ekappa = cosmology%ekappaofa(pert%a)
          pert%vis = pert%ekappa/pert%tauc
          pert%visdot = cosmology%vis%derivative(pert%a) * pert%a * pert%aH
          pert%Pdot = (pert%T2prime  - coop_sqrt6 *  pert%E2prime)/10.d0 * pert%aH
          pert%kchi = 1.d0 - pert%tau/cosmology%tau0
          pert%kchi = pert%k * cosmology%tau0 * (pert%kchi + exp(-1.d3*pert%kchi))
       endif

       pert%deltatr_mnu = 0.d0
       pert%deltap_mnu = 0.d0
    case(1)
       call coop_tbw("vector equations not written")
    case(2)
       O1_TEN_H_PRIME = O1_TEN_HPR
       if(pert%tight_coupling)then
          pert%T%F(2) =  (-16.d0/3.d0 )*O1_TEN_HPR*aHtauc
          pert%E%F(2) =  (-coop_sqrt6/4.d0)*pert%T%F(2) 
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
