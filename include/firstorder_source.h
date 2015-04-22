  subroutine coop_cosmology_firstorder_pert2source(this, pert, source, itau, ik)
    class(coop_cosmology_firstorder)::this
    type(coop_pert_object)::pert
    COOP_INT itau, ik
    type(coop_cosmology_firstorder_source)::source

    select case(source%m) 
    case(0) !!scalar
       source%saux(1, ik, itau) = (3.d0/8.d0)*pert%vis*pert%capP !*( (1.d0+ tanh(4.d0 - 20.d0*source%tau(itau)/this%tau0))/2.d0)  !!truncate at large z to get better numeric stability; this is not necessary for most smooth models.
       source%saux(2, ik, itau) = pert%O1_Phi+pert%O1_PSI
       source%saux(3, ik, itau) = pert%O1_PSI
       source%saux(6, ik, itau) = pert%O1_V_C
       source%saux(7, ik, itau) = pert%O1_V_B
       if(pert%has_rad_pert)then
          source%s(1,  ik, itau) =  (pert%O1_Phipr + pert%O1_PSIPR)*pert%aH*pert%ekappa    & !!ISW
               + pert%vis * (pert%O1_T(0)/4.d0 + pert%O1_Phi + pert%O1_V_B_PRIME/pert%kbyaH + pert%capP/8.0)  &
               + pert%visdot * (pert%O1_V_B/pert%k)
       else
          source%s(1,  ik, itau) =  (pert%O1_Phipr + pert%O1_PSIPR)*pert%aH*pert%ekappa    & !!ISW
               + pert%vis * (pert%O1_V_B_PRIME/pert%kbyaH + pert%capP/8.0)  &
               + pert%visdot * (pert%O1_V_B/pert%k)
       endif
       source%s(2, ik, itau) =pert%vis * pert%capP * (3.d0/8.d0)/ pert%kchi **2
       if(source%nsrc .ge.3)source%s(3, ik, itau) = -(pert%O1_Phi+pert%O1_PSI)*max(1.d0-source%chi(itau)/this%distlss, 0.d0)/max(source%chi(itau), 1.d-3)
       if(source%nsrc.ge.4) source%s(4, ik, itau) = -pert%vis 
    case(1) !!vector
       call coop_tbw("vector source to be done")
    case(2) !!tensor
       source%saux(1, ik, itau) = (coop_sqrt6/16.d0)*pert%vis*pert%capP
       source%s(1, ik, itau) = (pert%vis*pert%capP/4.d0 - pert%aH*pert%ekappa*pert%O1_TEN_HPR)*(coop_sqrt6/4.d0)/pert%kchi**2  !!tensor T
       source%s(2, ik, itau) = pert%vis*pert%capP * ((coop_sqrt6*3.d0/8.d0)/pert%kchi**2 - coop_sqrt6/16.d0) + (coop_sqrt6/4.d0)* (pert%vis*pert%Pdot + pert%visdot*pert%capP)/(pert%k*pert%kchi) 

       source%s(3, ik, itau) = pert%vis * pert%capP * ((coop_sqrt6/4.d0)/pert%kchi) + (coop_sqrt6/8.d0)*(pert%vis*pert%Pdot + pert%visdot*pert%capP)/pert%k
    end select
  end subroutine coop_cosmology_firstorder_pert2source
