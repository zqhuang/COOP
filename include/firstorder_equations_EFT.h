subroutine coop_cosmology_firstorder_equations(n, lna, y, yp, cosmology, pert)
    COOP_INT n
    type(coop_cosmology_firstorder)::cosmology
    type(coop_pert_object)::pert
    COOP_REAL, parameter::eps = 1.d-7
    COOP_REAL lna, y(0:n-1), yp(0:n-1)
    COOP_INT i, l, iq
    COOP_REAL a, aniso,  ktauc, ktaucdot, ktaucdd, aniso_prime, aHtauc, aHtau, aHsq, uterm, vterm, ma, doptdlna, M2a2H2, anisobyM2, anisobyM2_prime, u_prime, alpha_H_pp, wp1_de, w_de_prime, wp1eff_de, kbyaHsq_prime, sth, sth_prime, cmupp, auxt, auxt_prime, auxs, auxs_prime, psipr_th1, psipr_th2
    COOP_REAL :: pa2pr_g, pa2pr_nu
    COOP_REAL::  asq, Hsq
    COOP_INT, parameter::i_psipp = 1, i_mupp = 2, i_phi = 3, i_phip = 4, i_mup = 5, i_mu = 6, i_const = 7 , eq_phi = 1, eq_phip = 2, eq_psipp = 3, eq_mupp = 4, eq_aux = 5
    COOP_REAL::   pidot_c, mu_sol

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
    wp1_de = O0_DE(cosmology)%wp1ofa(a)
    wp1eff_de = O0_DE(cosmology)%wp1effofa(a)    
    w_de_prime = O0_DE(cosmology)%dwda(a)*a
    pert%pa2_de = (wp1_de - 1.d0 ) * pert%rhoa2_de
    

    
    pert%M2 = cosmology%mpsq(a)
    
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
    if(cosmology%f_alpha_H%initialized)then
       alpha_H_pp = a**2*cosmology%f_alpha_H%derivative2(a)+pert%alpha_H_prime
    else
       alpha_H_pp = 0.d0
    endif

    pert%rhoa2_mnu = 0.d0
    pert%pa2_mnu = 0.d0

    pert%rhoa2_matter = pert%rhoa2_g + pert%rhoa2_b + pert%rhoa2_c + pert%rhoa2_nu +pert%rhoa2_mnu
    pert%pa2_matter = pert%pa2_g + pert%pa2_nu + pert%pa2_mnu


    
    pert%rhoa2_sum = pert%rhoa2_matter + pert%rhoa2_de 
    pert%pa2_sum = pert%pa2_matter + pert%pa2_de
    aHsq = (pert%rhoa2_sum + cosmology%Omega_k())/3.d0/pert%M2  !!a^2H^2
    M2a2H2 = pert%M2 * aHsq
    pert%u = - (pert%rhoa2_de + pert%pa2_de)/(2.d0*M2a2h2)

    Hsq = aHsq/asq

    pert%aH = sqrt(aHsq)

    pert%daHdtau = -(pert%rhoa2_sum+3.d0*pert%pa2_sum)/6.d0/pert%M2
    pert%HdotbyHsq = pert%daHdtau/aHsq - 1.d0
    pert%HddbyH3 = cosmology%HddbyH3(a)
    pert%HdotbyHsq_prime = pert%HddbyH3 - 2.d0*pert%HdotbyHsq**2    

    u_prime = - (pert%alpha_M+2.d0*pert%HdotbyHsq)*pert%u + pert%rhoa2_de/(2.d0*M2a2h2)*(3.d0*wp1_de*wp1eff_de - w_de_prime)
    
    pa2pr_nu = O0_NU(cosmology)%dpa2da(a)*a
    pa2pr_g = O0_RADIATION(cosmology)%dpa2da(a)*a

    pert%R = 0.75d0 * pert%rhoa2_b/pert%rhoa2_g
    pert%tauc = cosmology%taucofa(a)
    pert%taucdot = cosmology%dot_tauc(a)

    ktauc = pert%k * pert%tauc
    ktaucdot = pert%k * pert%taucdot
    pert%kbyaH  = pert%k/pert%aH
    pert%kbyaHsq = pert%kbyaH**2
    kbyaHsq_prime = -2.d0 * pert%kbyaHsq * ( 1.d0 + pert%HdotbyHsq )
    aHtauc = pert%aH * pert%tauc
    doptdlna = 1.d0/aHtauc
    pert%tau = cosmology%tauofa(a)
    aHtau = pert%aH*pert%tau
    pert%latedamp = cosmology%late_damp_factor(pert%k, pert%tau)
    ktaucdd = pert%k*(cosmology%dot_tauc(a*1.01)-cosmology%dot_tauc(a/1.01))/(0.02/pert%aH)					      
    select case(pert%m)
    case(0)
       O1_PSI_PRIME = O1_PSIPR
       O1_DELTA_C_PRIME = - O1_V_C * pert%kbyaH + 3.d0 * O1_PSI_PRIME
       O1_DELTA_B_PRIME = - O1_V_B * pert%kbyaH + 3.d0 * O1_PSI_PRIME
       O1_T_PRIME(0) = (- O1_T(1) * pert%kbyaH/3.d0 + 4.d0 * O1_PSI_PRIME)*pert%latedamp
       O1_NU_PRIME(0) = (- O1_NU(1) * pert%kbyaH/3.d0 + 4.d0 * O1_PSI_PRIME)*pert%latedamp
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
       aniso = 0.6d0/pert%ksq * aniso
       anisobyM2 = aniso/pert%M2

       pert%deMat = 0.d0

       !!equation  for Phi (eq. 111)
       pert%deMat(i_phi, eq_phi) = 1.d0+pert%alpha_H
       pert%deMat(i_const, eq_phi) = anisobyM2- (1.d0+pert%alpha_T)*O1_PSI
       pert%deMat(i_mup, eq_phi) = - pert%alpha_H
       pert%deMat(i_mu, eq_phi)= (-pert%alpha_T + pert%alpha_M + pert%HdotbyHsq*pert%alpha_H)



       !!equation for Phi' (eq. 111 taking derivative)
       
       if(pert%tight_coupling)then
          aniso_prime = (pert%pa2_nu *(pert%kbyaH * (cosmology%klms_by_2lm1(2, 0, 0) *   O1_NU(1) - cosmology%klms_by_2lp1(3, 0, 0) *  O1_NU(3) ) )  + pa2pr_nu*O1_NU(2) +  pert%pa2_g * ( (8.d0/9.d0)*(ktauc*((O1_T(0) - 0.4d0*pert%T%F(2))*pert%kbyaH + 4.d0*pert%slip/aHtauc) + ktaucdot/pert%aH*O1_T(1)) ) + pa2pr_g * pert%T%F(2))*(0.6d0/pert%ksq)
          pert%deMat(i_phi, eq_phip) = pert%pa2_g * (32.d0/9.d0)* ktauc*pert%kbyaH/pert%M2*(0.6d0/pert%ksq)
       else
          aniso_prime = (pert%pa2_nu *(pert%kbyaH * (cosmology%klms_by_2lm1(2, 0, 0) *   O1_NU(1) - cosmology%klms_by_2lp1(3, 0, 0) *  O1_NU(3) ) )  + pa2pr_nu*O1_NU(2) +  pert%pa2_g * ( pert%kbyaH * (cosmology%klms_by_2lm1(2, 0, 0)*O1_T(1) - cosmology%klms_by_2lp1(3, 0, 0)*O1_T(3))  - (O1_T(2) - ((O1_T(2) - coop_sqrt6 * O1_E(2))/10.d0))/aHtauc) + pa2pr_g * pert%T%F(2))*(0.6d0/pert%ksq)
       endif
       anisobyM2_prime = aniso_prime/pert%M2 - anisobyM2*pert%alpha_M

       
       pert%deMat(i_phi, eq_phip) = pert%alpha_H_prime + pert%deMat(i_phi, eq_phip)
       pert%deMat(i_phip, eq_phip) = pert%deMat(i_phi, eq_phi)
       pert%deMat(i_mu,  eq_phip) = pert%alpha_M_prime - pert%alpha_T_prime + pert%HdotbyHsq*pert%alpha_H_prime + pert%HdotbyHsq_prime*pert%alpha_H
       pert%deMat(i_mup, eq_phip) = pert%deMat(i_mu, eq_phi) - pert%alpha_H_prime
       pert%deMat(i_mupp, eq_phip) = pert%deMat(i_mup, eq_phi)
       pert%deMat(i_const, eq_phip) =  anisobyM2_prime - (1.d0+pert%alpha_T)*O1_PSIPR - pert%alpha_T_prime * O1_PSI 
       


       !!equation for Psi'' (eq. 112)
       pert%deMat(i_psipp, eq_psipp) = 1.d0
       pert%deMat(i_mupp, eq_psipp) =  - pert%alpha_B
       pert%deMat(i_mup, eq_psipp) = pert%u - pert%alpha_B_prime - (4.d0+pert%alpha_M ) * pert%alpha_B + pert%alpha_K/6.d0 
       pert%deMat(i_mu, eq_psipp) =  (1.d0+pert%alpha_B) * pert%HdotbyHsq_prime &
            + pert%HdotbyHsq * (pert%alpha_B_prime+pert%HdotbyHsq - pert%u + 2.d0*pert%alpha_B - pert%alpha_K/6.d0 + (3.d0+pert%alpha_M+pert%HdotbyHsq)*(1.d0+pert%alpha_B))  &
            + pert%u - 2.d0*pert%pa2_matter/M2a2H2   &    !!here I am replacing p' with - 4p assuming no massive neutrinos
            + pert%kbyaHsq*(pert%alpha_H - pert%alpha_B)/3.d0
      pert%deMat(i_phi, eq_psipp) = 1.d0+2.d0*pert%HdotbyHsq - pert%u - pert%alpha_K/6.d0 + (2.d0+pert%HdotbyHsq)*pert%alpha_B + pert%alpha_B_prime + (3.d0+pert%alpha_M)*(1.d0+pert%alpha_B)
      pert%deMat(i_phip, eq_psipp) = 1.d0 + pert%alpha_B
      pert%deMat(i_const, eq_psipp)  =  - 0.5d0 * (O1_DELTA_B*pert%rhoa2_b/M2a2H2*(pert%cs2b - 1.d0/3.d0) + O1_DELTA_C * pert%rhoa2_c/M2a2H2*(-1.d0/3.d0)- 2.d0/3.d0 * pert%kbyaHsq * anisobyM2) &
               + (4.d0 + pert%HdotbyHsq + pert%alpha_M + pert%alpha_B)*O1_PSIPR &
               + ((1.d0 + pert%alpha_H)/3.d0*pert%kbyaHsq)*O1_PSI

      
      !!equations for mu''
       pert%deMat(i_psipp, eq_mupp) = 6.d0*pert%alpha_B
       pert%deMat(i_mupp, eq_mupp) = pert%alpha_K
       pert%deMat(i_mup, eq_mupp) = (3.d0+pert%alpha_M + pert%HdotbyHsq)*pert%alpha_K + pert%alpha_K_prime
       pert%deMat(i_mu, eq_mupp) = -pert%alpha_K * pert%HdotbyHsq_prime - pert%HdotbyHsq * ((3.d0 + pert%alpha_M + 2.d0*pert%HdotbyHsq)*pert%alpha_K + pert%alpha_K_prime ) &
            + 6.d0* (pert%HdotbyHsq*(pert%u + pert%alpha_B * (3.d0 + pert%alpha_M  + pert%HdotbyHsq*3.d0) + pert%alpha_B_prime) + pert%alpha_B*pert%HdotbyHsq_prime) &
            - 2.d0 * pert%kbyaHsq * ( pert%u + (pert%alpha_B - pert%alpha_H)*(1.d0+pert%alpha_M + pert%HdotbyHsq) +pert%alpha_T -pert%alpha_M +  pert%alpha_B_prime - pert%alpha_H_prime)
       pert%deMat(i_phip, eq_mupp) =  6.d0*pert%alpha_B- pert%alpha_K
       pert%deMat(i_phi, eq_mupp) = 6.d0*pert%u + (6.d0*pert%alpha_B - pert%alpha_K) * (3.d0 + pert%alpha_M) + 2.d0*(pert%alpha_B * 9.d0 - pert%alpha_K) * pert%HdotbyHsq  + 6.d0*pert%alpha_B_prime - pert%alpha_K_prime + 2.d0*pert%kbyaHsq*(pert%alpha_H - pert%alpha_B)
       pert%deMat(i_const, eq_mupp) =   (6.d0*(pert%alpha_B*(2.d0*pert%HdotbyHsq+3.d0+pert%alpha_M) + pert%u + pert%alpha_B_prime)+2.d0*pert%kbyaHsq*pert%alpha_H)*O1_PSIPR &
            + 2.d0*pert%kbyaHsq * (pert%alpha_M*(1.d0+pert%alpha_H) + pert%alpha_H  - pert%alpha_T + pert%alpha_H_prime)*O1_PSI

       !!Ok now try solving the equations
       pert%deMat(:, eq_phi) = pert%deMat(:, eq_phi)/pert%deMat(i_phi, eq_phi)
       pert%deMat(:, eq_phip) = pert%deMat(:, eq_phip)/pert%deMat(i_phip, eq_phip)
       pert%deMat(:, eq_psipp) = pert%deMat(:, eq_psipp)/pert%deMat(i_psipp, eq_psipp)
       pert%deMat(:, eq_phip) = pert%deMat(:, eq_phip) - pert%deMat(i_phi, eq_phip)*pert%deMat(:, eq_phi)
       pert%deMat(:, eq_psipp) = pert%deMat(:, eq_psipp) - pert%deMat(:, eq_phi)*pert%deMat(i_phi, eq_psipp) - pert%deMat(:, eq_phip) * pert%deMat(i_phip, eq_psipp)
       pert%deMat(:, eq_mupp) = pert%deMat(:, eq_mupp) - pert%deMat(:, eq_phi) * pert%deMat(i_phi, eq_mupp) - pert%deMat(:, eq_phip) * pert%deMat(i_phip, eq_mupp)
       pert%deMat(:, eq_mupp) = pert%deMat(:, eq_mupp) - pert%deMat(:, eq_psipp) * pert%deMat(i_psipp, eq_mupp)

       select case(pert%de_scheme)
       case(0) !!LCDM, in this case we don't really care about HPI
          O1_PHI = - pert%deMat(i_const, eq_phi)
          O1_PHI_PRIME = - pert%deMat(i_const, eq_phip) 
          O1_PSIPR_PRIME = - pert%deMat(i_const, eq_psipp) 
          O1_DE_HPI_PRIME = 0.d0
          O1_DE_HPIPR_PRIME = 0.d0
       case(1)
          if(abs(pert%deMat(i_mu, eq_mupp)) .lt. eps)then
             pert%deMat(i_mu, eq_mupp) = sign(eps,  pert%deMat(i_mu, eq_mupp))
          endif
          O1_DE_HPI = - pert%deMat(i_const, eq_mupp)/pert%deMat(i_mu, eq_mupp)
          pert%deMat(i_mup, eq_aux) =pert%deMat(i_mu, eq_mupp)
          pert%deMat(i_mu, eq_aux) = 6.d0 * ((pert%HdotbyHsq_prime+pert%alpha_T_prime - pert%alpha_M_prime)*pert%u +(pert%HdotbyHsq + pert%alpha_T - pert%alpha_M) * u_prime) - 2.d0*(pert%kbyaHsq*(u_prime + pert%alpha_T_prime - pert%alpha_M_prime) + kbyahsq_prime*(pert%u + pert%alpha_T - pert%alpha_M))
          pert%deMat(i_psipp, eq_aux) = 6.d0*pert%u
          pert%deMat(i_const, eq_aux) = (6.d0*u_prime+2.d0*pert%kbyaHsq*(pert%alpha_M - pert%alpha_T) + 6.d0*pert%u*(1.d0+pert%alpha_T) )*O1_PSIPR &
               + (2.d0*kbyahsq_prime*(pert%alpha_M - pert%alpha_T) + 2.d0*pert%kbyaHsq*(pert%alpha_M_prime - pert%alpha_T_prime) + 6.d0*(pert%alpha_T_prime*pert%u+(1.d0+pert%alpha_T)*u_prime))*O1_PSI &
               - 6.d0*(pert%u * anisobyM2_prime + u_prime * anisobyM2)
          pert%deMat(:, eq_psipp) =  pert%deMat(:, eq_psipp) - (pert%deMat(i_mup, eq_psipp)/pert%deMat(i_mup, eq_aux)) * pert%deMat(:, eq_aux)
          
          O1_PSIPR_PRIME = (- pert%deMat(i_const, eq_psipp)  - O1_DE_HPI *  pert%deMat(i_mu, eq_psipp))/pert%deMat(i_psipp, eq_psipp)
          O1_DE_HPIPR = -(pert%deMat(i_const, eq_aux)  + pert%deMat(i_mu, eq_aux)*O1_DE_HPI + pert%deMat(i_psipp, eq_aux)* O1_PSIPR_PRIME)/pert%deMat(i_mup, eq_aux)
          O1_PHI = - pert%deMat(i_const, eq_phi) - O1_DE_HPI * pert%deMat(i_mu, eq_phi) - O1_DE_HPIPR * pert%deMat(i_mup, eq_phi)
          O1_PHI_PRIME = - pert%deMat(i_const, eq_phip) - O1_DE_HPI * pert%deMat(i_mu, eq_phip) - O1_DE_HPIPR*pert%deMat(i_mup, eq_phip)

          O1_DE_HPI_PRIME = O1_DE_HPIPR
          O1_DE_HPIPR_PRIME = 0.d0
       case(2)
          if(abs(pert%deMat(i_mup, eq_mupp)) .lt. eps)then
             pert%deMat(i_mup, eq_mupp) = sign(eps,  pert%deMat(i_mup, eq_mupp))
          endif
          auxt = pert%alpha_M+pert%alpha_H *(1.d0+pert%alpha_M) - pert%alpha_T + pert%alpha_H_prime
          auxt_prime = pert%alpha_M_prime*(1.d0+pert%alpha_H) + pert%alpha_H_prime*(1.d0+pert%alpha_M)-pert%alpha_T_prime + alpha_H_pp
          auxs = 6.d0*pert%u + 2.d0*pert%kbyahsq*pert%alpha_H
          auxs_prime = 6.d0*u_prime + 2.d0*(kbyahsq_prime*pert%alpha_H + pert%kbyahsq*pert%alpha_H_prime)
          sth = auxs/(1.d0+pert%alpha_H)
          sth_prime = (auxs_prime - pert%alpha_H_prime/(1.d0+pert%alpha_H)*auxs)/(1.d0+pert%alpha_H)
          pert%deMat(i_mupp, eq_aux) = pert%deMat(i_mup, eq_mupp)
          pert%deMat(i_mup, eq_aux) = pert%deMat(i_mu, eq_mupp) &
               + pert%alpha_H_prime * sth + sth_prime*pert%alpha_H
          pert%deMat(i_mu, eq_aux) = 6.d0*(pert%HdotbyHsq_prime * pert%u + pert%HdotbyHsq*u_prime) &
               - 2.d0*kbyahsq_prime*(pert%u + pert%alpha_T - pert%alpha_M - pert%alpha_H_prime - pert%alpha_H*(1.d0+pert%alpha_M+pert%HdotbyHsq)) &
               - 2.d0*pert%kbyahsq*(u_prime + pert%alpha_T_prime - pert%alpha_M_prime - alpha_H_pp - pert%alpha_H_prime*(1.d0+pert%alpha_M+pert%HdotbyHsq) - pert%alpha_H*(pert%alpha_M_prime + pert%hdotbyhsq_prime)) &              
               - sth_prime * (pert%alpha_M - pert%alpha_T + pert%HdotbyHsq*pert%alpha_H) - (pert%alpha_M_prime - pert%alpha_T_prime + pert%HdotbyHsq_prime * pert%alpha_H + pert%alpha_H_prime*pert%HdotbyHsq)*sth
          pert%deMat(i_psipp, eq_aux)  = auxs
          pert%deMat(i_const, eq_aux) = ( 6.d0*u_prime+2.d0*(kbyahsq_prime*pert%alpha_H + pert%kbyahsq*pert%alpha_H_prime) + 2.d0*pert%kbyahsq*auxt + (1.d0+pert%alpha_T)*sth ) * O1_PSIPR &
               + (2.d0 * ( kbyahsq_prime*auxt + pert%kbyahsq*auxt_prime) + sth*pert%alpha_T_prime + sth_prime*(1.d0+pert%alpha_T))*O1_PSI &
               - sth*anisobyM2_prime - sth_prime*anisobyM2

          pert%deMat(:, eq_psipp) = pert%deMat(:, eq_psipp) - (pert%deMat(i_mupp, eq_psipp)/pert%deMat(i_mupp, eq_aux))*pert%deMat(:, eq_aux)
          O1_DE_HPIPR = (-pert%deMat(i_const, eq_psipp)-pert%deMat(i_mu, eq_psipp)*O1_DE_HPI)/pert%deMat(i_mup, eq_psipp)
          O1_DE_HPI_PRIME = O1_DE_HPIPR
          O1_PHI = - pert%deMat(i_const, eq_phi) - pert%deMat(i_mu,eq_phi)*O1_DE_HPI - pert%deMat(i_mup, eq_phi) * O1_DE_HPIPR
          psipr_th1 = (- O1_PHI - O1_PSI*(2.d0*pert%kbyahsq/auxs*auxt) - O1_DE_HPI*(pert%HdotbyHsq+2.d0*pert%kbyahsq/auxs*(auxt-pert%u)))
          psipr_th2 = (- pert%rhoa2_b*O1_DELTA_B &
               - pert%rhoa2_c*O1_DELTA_C &
               - pert%rhoa2_g*O1_T(0) &
               - pert%rhoa2_nu*O1_NU(0))/m2a2h2/6.d0 &             
               - O1_PHI-pert%kbyahsq/3.d0*(1.d0+pert%alpha_H)*O1_PSI - O1_DE_HPI*(pert%u+pert%kbyahsq*pert%alpha_H/3.d0)
          O1_PSIPR_PRIME = (- O1_PSIPR + (psipr_th1+psipr_th2*pert%kbyahsq)/(1.d0+pert%kbyahsq))*1.d4
          O1_DE_HPIPR_PRIME = 0.d0
          O1_PHI_PRIME = 0.d0
       case(3)
          if(abs(pert%deMat(i_mupp, eq_mupp)) .lt. eps)then
             pert%deMat(i_mupp, eq_mupp) = sign(eps,  pert%deMat(i_mupp, eq_mupp))
          endif
          O1_DE_HPI_PRIME = O1_DE_HPIPR          
          O1_DE_HPIPR_PRIME = (-pert%deMat(i_const, eq_mupp) - pert%deMat(i_mu, eq_mupp)*O1_DE_HPI - pert%deMat(i_mup, eq_mupp)*O1_DE_HPIPR)/pert%deMat(i_mupp, eq_mupp)
          O1_PSIPR_PRIME = - pert%deMat(i_const, eq_psipp) - pert%deMat(i_mu, eq_psipp)*O1_DE_HPI - pert%deMat(i_mup, eq_psipp)*O1_DE_HPIPR - pert%deMat(i_mupp, eq_psipp)*O1_DE_HPIPR_PRIME
          O1_PHI = - pert%deMat(i_const, eq_phi) - pert%deMat(i_mu,eq_phi)*O1_DE_HPI - pert%deMat(i_mup, eq_phi) * O1_DE_HPIPR
          O1_PHI_PRIME = -pert%deMat(i_const, eq_phip) - pert%deMat(i_mu, eq_phip)*O1_DE_HPI - pert%deMat(i_mup,eq_phip)*O1_DE_HPIPR - pert%deMat(i_mupp, eq_phip)*O1_DE_HPIPR_PRIME
       case default
          stop "unknown EFT DE scheme"
       end select
       pert%delta_gamma = pert%latedamp * O1_T(0) - 4.d0*(1.d0-pert%latedamp)*(O1_PHI + O1_V_B/ktauc)
       
       !!velocities
       O1_V_C_PRIME = - O1_V_C + pert%kbyaH * O1_PHI
       O1_NU_PRIME(1) = ((O1_NU(0) + 4.d0*O1_PHI - 0.4d0 * O1_NU(2))*pert%kbyaH)*pert%latedamp
       O1_V_B_PRIME = - O1_V_B + pert%kbyaH * (O1_PHI + pert%cs2b * O1_DELTA_B) - pert%slip/(pert%R * aHtauc)
       O1_T_PRIME(1) = ((O1_T(0) + 4.d0*O1_PHI - 0.4d0*pert%T%F(2))*pert%kbyaH + 4.d0*pert%slip/aHtauc)*pert%latedamp


       !!higher moments
       !!massless neutrinos
       do l = 2, pert%nu%lmax - 1
          O1_NU_PRIME(l) =  (pert%kbyaH * (cosmology%klms_by_2lm1(l, 0, 0) *   O1_NU( l-1 ) - cosmology%klms_by_2lp1(l+1, 0, 0) *  O1_NU( l+1 ) ))*pert%latedamp
       enddo

       O1_NU_PRIME(pert%nu%lmax) = (pert%kbyaH * (pert%nu%lmax +0.5d0)/(pert%nu%lmax-0.5d0)*  O1_NU(pert%nu%lmax-1) -  (pert%nu%lmax+1)* O1_NU(pert%nu%lmax)/(aHtau))*pert%latedamp

       if(pert%tight_coupling)then

          pert%T2prime = (8.d0/9.d0)*(ktauc*O1_T_PRIME(1) + ktaucdot/pert%aH*O1_T(1)) !2.d0*pert%T%F(2)
          pert%E%F(2) = -coop_sqrt6/4.d0 * pert%T%F(2)
          pert%E2prime = -coop_sqrt6/4.d0 * pert%T2prime 
          pert%capP = (pert%T%F(2) - coop_sqrt6 * pert%E%F(2))/10.d0

       else
          !!T
          pert%capP = (O1_T(2) - coop_sqrt6 * O1_E(2))/10.d0
          O1_T_PRIME(2) =  (pert%kbyaH * (cosmology%klms_by_2lm1(2, 0, 0)*O1_T(1) - cosmology%klms_by_2lp1(3, 0, 0)*O1_T(3))  - (O1_T(2) - pert%capP)/aHtauc)*pert%latedamp
          pert%T2prime =  O1_T_PRIME(2)
          do l = 3, pert%T%lmax -1
             O1_T_PRIME(l) = (pert%kbyaH * (cosmology%klms_by_2lm1(l, 0, 0)*O1_T(l-1) -cosmology%klms_by_2lp1(l+1, 0, 0)*O1_T(l+1))  - O1_T(l)/aHtauc)*pert%latedamp
          enddo
          O1_T_PRIME(pert%T%lmax) =  (pert%kbyaH *((pert%T%lmax+0.5d0)/(pert%T%lmax-0.5d0))*  O1_T(pert%T%lmax-1) &
               -  ((pert%T%lmax+1)/aHtau + doptdlna) * O1_T(pert%T%lmax))*pert%latedamp
          !!E
          O1_E_PRIME(2) = (pert%kbyaH * ( - cosmology%klms_by_2lp1(3, 0, 2)*O1_E(3))  - (O1_E(2) + coop_sqrt6 * pert%capP)/aHtauc)*pert%latedamp
          pert%E2prime = O1_E_PRIME(2)

          do l = 3, pert%E%lmax - 1
             O1_E_PRIME(l) = (pert%kbyaH * (cosmology%klms_by_2lm1(l, 0, 2)*O1_E(l-1) - cosmology%klms_by_2lp1(l+1, 0, 2)*O1_E(l+1)) - O1_E(l)/aHtauc)*pert%latedamp
          enddo
          O1_E_PRIME(pert%E%lmax) =  (pert%kbyaH *( (pert%E%lmax-2.d0/pert%E%lmax+0.5d0)/(pert%E%lmax-2.d0/pert%E%lmax-0.5d0)) *  O1_E(pert%E%lmax-1) &
               -  ((pert%E%lmax-2.d0/pert%E%lmax+1)/aHtau + doptdlna) * O1_E(pert%E%lmax))*pert%latedamp

       endif
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
          pert%T%F(2) =O1_T(2)
          pert%E%F(2) = O1_E(2)
       endif
       pert%capP = (pert%T%F(2) - coop_sqrt6 * pert%E%F(2))/10.d0
       aniso = pert%pa2_g * pert%T%F(2) + pert%pa2_nu * O1_NU(2)  
       O1_TEN_HPR_PRIME = -(2.d0+pert%daHdtau/aHsq)*O1_TEN_HPR &
            - pert%kbyaH**2 * O1_TEN_H &
            + 0.4d0/aHsq * aniso/pert%M2

       O1_NU_PRIME(2) =  (pert%kbyaH * ( - cosmology%klms_by_2lp1(3, 2, 0) *  O1_NU( 3 ) ) - 4.d0*O1_TEN_HPR)*pert%latedamp
       do l = 3, pert%nu%lmax - 1
          O1_NU_PRIME(l) =  (pert%kbyaH * (cosmology%klms_by_2lm1(l, 2, 0) *   O1_NU( l-1 ) - cosmology%klms_by_2lp1(l+1, 2, 0) *  O1_NU( l+1 ) ))*pert%latedamp
       enddo

       O1_NU_PRIME(pert%nu%lmax) = (pert%kbyaH * (pert%nu%lmax-2.d0/pert%nu%lmax+0.5d0)/(pert%nu%lmax-2.d0/pert%nu%lmax-0.5d0)*  O1_NU(pert%nu%lmax-1) &
            -  (pert%nu%lmax-2.d0/pert%nu%lmax+1)* O1_NU(pert%nu%lmax)/(aHtau))*pert%latedamp

       if(pert%tight_coupling)then
          pert%T2prime = -16.d0/3.d0*(O1_TEN_HPR_PRIME*aHtauc - O1_TEN_HPR * (pert%taucdot  + pert%tauc*pert%daHdtau/pert%aH)/aHtauc**2)
          pert%E2prime = (-coop_sqrt6/4.d0)*pert%T2prime
       else
          O1_T_PRIME(2) =  (pert%kbyah * ( -  cosmology%klms_by_2lp1(3, 2, 0) * O1_T(3))  -  ( O1_T(2) - pert%capP)/aHtauc &
               - O1_TEN_HPR*4.d0)*pert%latedamp
          pert%T2prime = O1_T_PRIME(2)
          O1_E_PRIME(2) = (pert%kbyah * ( - cosmology%fourbyllp1(2)*O1_B(2)- cosmology%klms_by_2lp1(3, 2, 2) * O1_E(3)) - (O1_E(2) + coop_sqrt6 * pert%capP)/aHtauc)*pert%latedamp
          pert%E2prime = O1_E_PRIME(2)
          O1_B_PRIME(2) = (pert%kbyah * ( cosmology%fourbyllp1(2)*O1_E(2)- cosmology%klms_by_2lp1(3, 2, 2) * O1_B(3)) - O1_B(2)/aHtauc)*pert%latedamp
          do l = 3, pert%T%lmax -1 
             O1_T_PRIME(l) = (pert%kbyah * (cosmology%klms_by_2lm1(l, 2, 0) * O1_T(l-1) -  cosmology%klms_by_2lp1(l+1, 2, 0) * O1_T(l+1)) -  O1_T(l)/aHtauc)*pert%latedamp
          enddo
          O1_T_PRIME(pert%T%lmax) = (pert%kbyah  *((pert%T%lmax-2.d0/pert%T%lmax+0.5d0)/(pert%T%lmax-2.d0/pert%T%lmax-0.5d0)) * O1_T(pert%T%lmax-1) &
               - (doptdlna + (pert%T%lmax-2.d0/pert%T%lmax+1)/aHtau)* O1_T(pert%T%lmax))*pert%latedamp
          if(pert%E%lmax .ne. pert%B%lmax) stop "lmax(E) must be the same as lmax(B)"
          do l=3, pert%E%lmax -1
             O1_E_PRIME(l) = (pert%kbyah * (cosmology%klms_by_2lm1(l, 2, 2) * O1_E(l-1) &
                  - cosmology%fourbyllp1(l)*O1_B(l) &
                  - cosmology%klms_by_2lp1(l+1, 2, 2) * O1_E(l+1)) - O1_E(l)/aHtauc)*pert%latedamp
             O1_B_PRIME(l) = (pert%kbyah * (cosmology%klms_by_2lm1(l, 2, 2) * O1_B(l-1) + cosmology%fourbyllp1(l)*O1_E(l)- cosmology%klms_by_2lp1(l+1, 2, 2) * O1_B(l+1)) -  O1_B(l)/aHtauc)*pert%latedamp
          enddo
          O1_E_PRIME(pert%E%lmax) = (pert%kbyah  *((pert%E%lmax-4.d0/pert%E%lmax+0.5d0)/(pert%E%lmax-4.d0/pert%E%lmax-0.5d0)) * O1_E(pert%E%lmax-1) &
               - cosmology%fourbyllp1(pert%E%lmax)*O1_B(pert%E%lmax) &
               - (doptdlna + (pert%E%lmax-4.d0/pert%E%lmax+1)/aHtau)* O1_E(pert%E%lmax))*pert%latedamp
          O1_B_PRIME(pert%E%lmax) = (pert%kbyah  *((pert%E%lmax-4.d0/pert%E%lmax+0.5d0)/(pert%E%lmax-4.d0/pert%E%lmax-0.5d0)) * O1_B(pert%E%lmax-1)  + cosmology%fourbyllp1(pert%E%lmax)*O1_E(pert%E%lmax) &
               - (doptdlna + (pert%E%lmax-4.d0/pert%E%lmax+1)/aHtau)* O1_B(pert%E%lmax))*pert%latedamp
       endif
       if(pert%want_source)then
          pert%ekappa = cosmology%ekappaofa(pert%a)
          pert%vis = pert%ekappa/pert%tauc
          pert%visdot = cosmology%vis%derivative(pert%a) * pert%a * pert%aH
          pert%Pdot = (pert%T2prime  - coop_sqrt6 *  pert%E2prime)/10.d0 * pert%aH
          pert%kchi = 1.d0 - pert%tau/cosmology%tau0
          pert%kchi = pert%k * cosmology%tau0 * (pert%kchi + exp(-1.d3*pert%kchi))
       endif
    case default
       call coop_return_error("firstorder_equations", "Unknown m = "//trim(coop_num2str(pert%m)), "stop")
    end select
    pert%capP = pert%capP*pert%latedamp    
  end subroutine coop_cosmology_firstorder_equations
