  subroutine coop_cosmology_firstorder_pert2source(this, pert, source, itau, ik)
    class(coop_cosmology_firstorder)::this
    type(coop_pert_object)::pert
    COOP_INT itau, ik
    type(coop_cosmology_firstorder_source)::source

    select case(source%m) 
    case(0) !!scalar
       source%saux(1, ik, itau) = (3.d0/8.d0)*pert%vis*pert%capP
       source%saux(coop_aux_index_Weyl, ik, itau) = pert%O1_Phi+pert%O1_PSI
       source%saux(coop_aux_index_psi, ik, itau) = pert%O1_PSI
#if DO_COUPLED_DE
    source%saux(coop_aux_index_delta_sync, ik, itau) = (pert%O1_DELTA_C - pert%O1_V_C/pert%kbyaH * O0_CDM(this)%dlnrhodlna(pert%a))/pert%k**2        
#else       
    source%saux(coop_aux_index_delta_sync, ik, itau) = (pert%O1_DELTA_C + 3.d0*pert%O1_V_C/pert%kbyaH)/pert%k**2
#endif       
       source%s(coop_index_source_T,  ik, itau) =  (pert%O1_Phipr + pert%O1_PSIPR)*pert%aH*pert%ekappa    & !!ISW
            + pert%vis * (pert%delta_gamma/4.d0 + pert%O1_Phi + pert%O1_V_B_PRIME/pert%kbyaH + pert%capP/8.0)  &
            + pert%visdot * (pert%O1_V_B/pert%k)
       
       source%s(coop_index_source_E, ik, itau) =pert%vis * pert%capP * (3.d0/8.d0)/ pert%kchi **2
       source%s(coop_index_source_Len, ik, itau) = -(pert%O1_Phi+pert%O1_PSI)*max(1.d0-source%chi(itau)/this%distlss, 0.d0)/max(source%chi(itau), 0.005*source%distlss)
!    source%s(coop_index_source_B, ik, itau) = 0.
#if DO_ZETA_TRANS
       if(coop_zeta_user_specified_weight%initialized)then
          source%s(coop_index_source_zeta, ik, itau) = -coop_zeta_user_specified_weight%eval(source%chi(itau))
       elseif(coop_zeta_single_slice_chi .gt. 0.d0)then
          source%s(coop_index_source_zeta, ik, itau) = 0.d0          
       else
          source%s(coop_index_source_zeta, ik, itau) = pert%vis
       endif
#endif       
    case(1) !!vector
       call coop_tbw("vector source to be done")
    case(2) !!tensor
       source%saux(1, ik, itau) = (coop_sqrt6/16.d0)*pert%vis*pert%capP
       source%s(coop_index_source_T, ik, itau) = (pert%vis*pert%capP/4.d0 - pert%aH*pert%ekappa*pert%O1_TEN_HPR)*(coop_sqrt6/4.d0)/pert%kchi**2  !!tensor T
       source%s(coop_index_source_E, ik, itau) = pert%vis*pert%capP * ((coop_sqrt6*3.d0/8.d0)/pert%kchi**2 - coop_sqrt6/16.d0) + (coop_sqrt6/4.d0)* (pert%vis*pert%Pdot + pert%visdot*pert%capP)/(pert%k*pert%kchi) 

       source%s(coop_index_source_B, ik, itau) = pert%vis * pert%capP * ((coop_sqrt6/4.d0)/pert%kchi) + (coop_sqrt6/8.d0)*(pert%vis*pert%Pdot + pert%visdot*pert%capP)/pert%k
    end select
  end subroutine coop_cosmology_firstorder_pert2source
