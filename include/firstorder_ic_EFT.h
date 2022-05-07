  subroutine coop_cosmology_firstorder_set_initial_conditions(this, pert, m, k, tau)
    class(coop_cosmology_firstorder)::this
    class(coop_pert_object)::pert
    COOP_INT :: m 
    COOP_REAL tau, k, Rnu, M2
    COOP_REAL:: mucoef, phicoef, a

    pert%k = k
    pert%tight_coupling = .true.
    call pert%init(m = m, nu_mass = this%mnu_by_Tnu, de_genre = this%de_genre)
    call pert%set_zero()
    
    if(this%index_massivenu .ne. 0)then
       pert%num_mnu_ratio = O0_MASSIVENU(this)%Omega_massless/O0_NU(this)%Omega
    else
       pert%num_mnu_ratio = 0.d0
    endif

    select case(trim(pert%initial_conditions))
    case("adiabatic")
       select case(pert%m)
       case(0)
          Rnu = this%Omega_nu / this%Omega_r
          pert%O1_Phi = -coop_primordial_zeta_norm/(1.5d0 + 0.4d0*Rnu)
          pert%O1_PhiPr = 0.d0
          pert%O1_PSI =  (1.+0.4*Rnu)*pert%O1_Phi
          pert%O1_PSIPR = 0.d0

          pert%O1_DELTA_B = -1.5d0*pert%O1_Phi
          pert%O1_DELTA_C = pert%O1_DELTA_B
          pert%O1_T(0) =  -2.d0*pert%O1_Phi
          pert%O1_NU(0) =  pert%O1_T(0)
          

          pert%O1_V_B = pert%O1_Phi/2.d0*k*tau
          pert%O1_V_C = pert%O1_V_B
          pert%O1_NU(1) = 4.d0*pert%O1_V_C
          pert%O1_T(1) = pert%O1_NU(1)
          
          pert%O1_NU(2) = (2.d0/3.d0) * pert%O1_Phi * (k*tau)**2

          a = this%aoftau(tau)
          pert%alpha_M = this%alpha_M(a)
          pert%alpha_B = this%alpha_B(a)
          pert%alpha_K = this%alpha_K(a)
          pert%alpha_H = this%alpha_H(a)

          pert%alpha_M_prime = this%alpha_M_prime(a)
          pert%alpha_B_prime = this%alpha_B_prime(a)
          pert%alpha_K_prime = this%alpha_K_prime(a)
          pert%alpha_H_prime = this%alpha_H_prime(a)

          mucoef = 2.d0 * ((pert%alpha_M - 1.d0)*pert%alpha_K + pert%alpha_K_prime ) &
               + 6.d0* (-2.d0 *(pert%alpha_B * (-3.d0 + pert%alpha_M) + pert%alpha_B_prime) ) 


          if(abs(mucoef) .gt. 1.d-30)then
             phicoef  = (6.d0*pert%alpha_B - pert%alpha_K) * (3.d0 + pert%alpha_M) -4.d0 *(pert%alpha_B * 9.d0 - pert%alpha_K)  + 6.d0*pert%alpha_B_prime - pert%alpha_K_prime              
             pert%O1_DE_HPI = (-phicoef/mucoef)*pert%O1_Phi
          else
             pert%O1_DE_HPI = pert%O1_Phi/2.d0
          endif
          pert%O1_DE_HPIPR = 0.d0
       case(1)
          call coop_tbw("vector initialization")
       case(2)
          pert%O1_TEN_H = -coop_primordial_zeta_norm/coop_sqrt6
          pert%O1_TEN_HPR = -(k*tau)**2/5.d0 
       end select
    case default
       stop "unknown initial conditions"
    end select
  end subroutine coop_cosmology_firstorder_set_initial_conditions

