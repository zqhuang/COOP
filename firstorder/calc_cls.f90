program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

#if DO_ZETA_TRANS
  logical::zeta_single_slice = .false.
#endif  
  !!----------------------------------------
  COOP_STRING::params_file, output
  type(coop_dictionary)::params
  type(coop_cosmology_firstorder)::cosmology
  COOP_REAL::norm, lnorm
  COOP_INT::l
  type(coop_file)::fp
  logical success
  !!----------------------------------------
  if(iargc() .lt. 1)then
     write(*,*) "Syntax:"
     write(*,*) "./CalcCls  params.ini"
     stop
  endif
  coop_feedback_level = 4 
  call coop_get_Input(1, params_file)
  call coop_load_dictionary(params_file, params)
  call cosmology%init_from_dictionary(params, level = coop_init_level_set_tens, success = success)
  if(.not. success) stop "Unhealthy model: perturbations blow up exponentially."
  print*, "distance to last scattering surface: ", cosmology%distlss/cosmology%H0Mpc(), " Mpc"
  print*, "distance to z=1089: ", cosmology%comoving_distance(1.d0/1090.0)/cosmology%H0Mpc(), " Mpc"
  if(abs(cosmology%mpsq0-1).gt. 1.d-6) print*, "M*^2 = ", cosmology%mpsq0

  call coop_dictionary_lookup(params, "root", output, default_val="test")
  
  !!----------------------------------------
#if DO_ZETA_TRANS
  if(zeta_single_slice)coop_zeta_single_slice_chi = cosmology%distlss
#endif  
  norm =COOP_DEFAULT_TCMB**2*1.d12 !cosmology%Tcmb()**2*1.d12    
  call fp%open(trim(output)//"_scalCls.dat","w")
  write(*,*) "scalar Cls are saved  in "//trim(output)//"_scalCls.dat"  
  write(fp%unit, "(A8, 5A16)") "# ell ", "   TT  ",  "   EE  ",  "   TE   ", " PhiPhi ", " TPhi "

  do l=cosmology%source(0)%trans%lmin, cosmology%source(0)%trans%lmax
     lnorm =  l*(l+1.d0)/coop_2pi*norm
     write(fp%unit, "(I8, 6E16.7)") l, cosmology%source(0)%Cls(coop_index_ClTT, l)*lnorm,  cosmology%source(0)%Cls(coop_index_ClEE, l)*lnorm,  cosmology%source(0)%Cls(coop_index_ClTE, l)*lnorm, cosmology%source(0)%Cls(coop_index_ClLenLen, l)*norm*(l*(l+1.d0))**2, cosmology%source(0)%Cls(coop_index_ClTLen, l)*norm*(l*(l+1.d0))**1.5
  enddo
  call fp%close()
  call fp%open(trim(output)//"_lensedCls.dat","w")
  write(*,*) "lensed Cls are saved  in "//trim(output)//"_lensedCls.dat"     
  write(fp%unit, "(A8, 5A16)") "# ell ", "   TT  ",  "   EE  ",  " BB  ", "  TE "
  do l=cosmology%source(0)%trans%lmin, cosmology%source(0)%trans%lmax
     lnorm = l*(l+1.d0)/coop_2pi*norm
     write(fp%unit, "(I5, 20E16.7)") l, cosmology%source(0)%cls_lensed(coop_index_ClTT, l)*lnorm, cosmology%source(0)%cls_lensed(coop_index_ClEE, l)*lnorm,  cosmology%source(0)%cls_lensed(coop_index_ClBB, l)*lnorm,  cosmology%source(0)%cls_lensed(coop_index_ClTE, l)*lnorm
  enddo
  call fp%close()
  if(cosmology%has_tensor)then  !! do tensor modes
     write(*,*) "tensor Cls are saved  in "//trim(output)//"_tensCls.dat"     
    call fp%open(trim(output)//"_tensCls.dat","w")
     write(fp%unit, "(A8, 5A16)") "# ell ", "   TT  ",  "   EE  ",  " BB  ", "  TE "
     do l=cosmology%source(2)%trans%lmin, cosmology%source(2)%trans%lmax
        lnorm =  l*(l+1.d0)/coop_2pi*norm
        write(fp%unit, "(I8, 6E16.7)") l, cosmology%source(2)%cls(coop_index_ClTT, l)*lnorm,  cosmology%source(2)%cls(coop_index_ClEE, l)*lnorm,  cosmology%source(2)%cls(coop_index_ClBB, l)*lnorm,  cosmology%source(2)%cls(coop_index_ClTE, l)*lnorm
     enddo
     call fp%close()
     call fp%open(trim(output)//"_totCls.dat","w")
     write(*,*) "total Cls are saved  in "//trim(output)//"_totCls.dat"     
     write(fp%unit, "(A8, 5A16)") "# ell ", "   TT  ",  "   EE  ",  " BB  ", "  TE "
     do l=max(cosmology%source(2)%trans%lmin, cosmology%source(0)%trans%lmin), min(cosmology%source(0)%trans%lmax, cosmology%source(2)%trans%lmax)
        lnorm = l*(l+1.d0)/coop_2pi*norm
        write(fp%unit, "(I5, 20E16.7)") l, (cosmology%source(0)%cls_lensed(coop_index_ClTT, l)+cosmology%source(2)%cls_lensed(coop_index_ClTT, l))*lnorm,  &
             (cosmology%source(0)%cls_lensed(coop_index_ClEE, l)+cosmology%source(2)%cls_lensed(coop_index_ClEE, l))*lnorm,  &
             (cosmology%source(0)%cls_lensed(coop_index_ClBB, l)+cosmology%source(2)%cls_lensed(coop_index_ClBB, l))*lnorm, &
             (cosmology%source(0)%cls_lensed(coop_index_ClTE, l)+cosmology%source(2)%cls_lensed(coop_index_ClTE, l))*lnorm
     enddo
     call fp%close()
     
  endif

contains

  subroutine generate_latevis()
    COOP_INT, parameter::n = 1024
    COOP_REAL::chi(n), vis(n), a, chiend
    COOP_INT::i
#if DO_ZETA_TRANS    
    call coop_set_uniform(n, chi, 0.d0, cosmology%tau0)
    chiend = cosmology%tau0 - cosmology%dadtau(1.d0/(1.d0+cosmology%zre+cosmology%deltaz*5.d0))
    do i=1, n
       if(chi(i) .lt. chiend)then
          a = cosmology%aoftau(cosmology%tau0 - chi(i))
          vis(i) = cosmology%visofa(a)
       else
          vis(i) = 0.d0
       endif
    enddo
    vis = vis/(sum(vis)*(chi(2)-chi(1)))
    call coop_zeta_user_specified_weight%init(n = n, xmin = chi(1), xmax = chi(n), f = vis, method = COOP_INTERPOLATE_LINEAR, name="latevis")
#endif    
  end subroutine generate_latevis


  subroutine generate_earlyvis()
    COOP_INT, parameter::n = 1024
    COOP_REAL::chi(n), vis(n), a, chiend
    COOP_INT::i
#if DO_ZETA_TRANS    
    call coop_set_uniform(n, chi, 0.d0, cosmology%tau0)
    chiend = cosmology%tau0 - cosmology%tauofa(1.d0/(1.d0+cosmology%zre+cosmology%deltaz*5.d0))
    do i=1, n
       if(chi(i) .ge. chiend)then
          a = cosmology%aoftau(cosmology%tau0 - chi(i))
          vis(i) = cosmology%visofa(a)
       else
          vis(i) = 0.d0
       endif
    enddo
    vis = vis/(sum(vis)*(chi(2)-chi(1)))
    call coop_zeta_user_specified_weight%init(n = n, xmin = chi(1), xmax = chi(n), f = vis, method = COOP_INTERPOLATE_LINEAR, name="earlyvis")
#endif    
  end subroutine generate_earlyvis

end program test
