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
  COOP_INT, parameter::lmin = 2, lmax = 2608
  COOP_REAL::Cls(coop_num_Cls, lmin:lmax), tensCls(coop_num_Cls, lmin:lmax), lensedCls(coop_num_Cls, lmin:lmax), ells(lmin:lmax)
  COOP_REAL::norm, lnorm, lambda
  COOP_INT::l
  logical success
  type(coop_file)::fp
  !!----------------------------------------
  if(iargc() .lt. 1)then
     write(*,*) "Syntax:"
     write(*,*) "./CalcCls  params.ini"
     stop
  endif
  
  call coop_get_Input(1, params_file)
  call coop_load_dictionary(params_file, params)
  call cosmology%init_from_dictionary(params)
  call coop_dictionary_lookup(params, "root", output, default_val="test")
  
  !!----------------------------------------
#if DO_ZETA_TRANS
  if(zeta_single_slice)coop_zeta_single_slice_chi = cosmology%distlss
#endif  
  call cosmology%init_source(0)
  call cosmology%compute_source(0, success = success)
  if(.not. success) stop "Solution blows up exponentially; Model is ruled out."
  write(*,*) "linear perturbation done: sigma_8  = ", cosmology%sigma_8  
  call cosmology%source(0)%get_all_cls(lmin, lmax, Cls)
  norm =COOP_DEFAULT_TCMB**2*1.d12 !cosmology%Tcmb()**2*1.d12    
  call fp%open(trim(output)//"_scalCls.dat","w")
  write(*,*) "scalar Cls are saved  in "//trim(output)//"_scalCls.dat"  
  write(fp%unit, "(A8, 5A16)") "# ell ", "   TT  ",  "   EE  ",  "   TE   ", " PhiPhi ", " TPhi "
  do l = lmin, lmax
     lnorm =  l*(l+1.d0)/coop_2pi*norm
     write(fp%unit, "(I8, 6E16.7)") l, Cls(coop_index_ClTT, l)*lnorm,  Cls(coop_index_ClEE, l)*lnorm,  Cls(coop_index_ClTE, l)*lnorm, Cls(coop_index_ClLenLen, l)*norm*(l*(l+1.d0))**2, Cls(coop_index_ClTLen, l)*norm*(l*(l+1.d0))**1.5
  enddo
  call fp%close()
  call coop_get_lensing_Cls(lmin, lmax, Cls, lensedCls)
  lensedCls = lensedCls + Cls
  call fp%open(trim(output)//"_lensedCls.dat","w")
  write(*,*) "lensed Cls are saved  in "//trim(output)//"_lensedCls.dat"     
  write(fp%unit, "(A8, 5A16)") "# ell ", "   TT  ",  "   EE  ",  " BB  ", "  TE "
  do l=lmin, lmax
     lnorm = l*(l+1.d0)/coop_2pi*norm
     write(fp%unit, "(I5, 20E16.7)") l, lensedCls(coop_index_ClTT, l)*lnorm, lensedCls(coop_index_ClEE, l)*lnorm,  lensedCls(coop_index_ClBB, l)*lnorm,  lensedCls(coop_index_ClTE, l)*lnorm
  enddo
  call fp%close()
#if DO_ZETA_TRANS
  call fp%open(trim(output)//"_zetaCls.dat","w")
  write(*,*) "zeta Cls are saved  in "//trim(output)//"_zetaCls.dat"     
  write(fp%unit, "(A7, 5A16)") "# ell", "  TT  ",  " EE  ",  " ZZ ", " TE ", " TZ ",  " EZ "
  do l=lmin, lmax
     ells(l) = l
     lnorm = l*(l+1.d0)/coop_2pi*norm
     if(zeta_single_slice)then
        write(fp%unit, "(I5, 20E16.7)") l, Cls(coop_index_ClTT, l)*lnorm, Cls(coop_index_ClEE, l)*lnorm,  Cls(coop_index_Clzetazeta, l)*lnorm, Cls(coop_index_ClTE, l)*lnorm, Cls(coop_index_ClTzeta, l)*lnorm, Cls(coop_index_ClEzeta, l)*lnorm, cosmology%Clzetazeta_at_R(l, cosmology%distlss)*lnorm        
     else
        write(fp%unit, "(I5, 20E16.7)") l, Cls(coop_index_ClTT, l)*lnorm, Cls(coop_index_ClEE, l)*lnorm,  Cls(coop_index_Clzetazeta, l)*lnorm, Cls(coop_index_ClTE, l)*lnorm, Cls(coop_index_ClTzeta, l)*lnorm, Cls(coop_index_ClEzeta, l)*lnorm
     endif
  enddo
  call fp%close()
#endif  
  if(cosmology%r .gt. 0.d0)then  !! do tensor modes
     call cosmology%compute_source(2, success = success)
     call cosmology%source(2)%get_all_cls(lmin, lmax, tensCls)
     write(*,*) "tensor Cls are saved  in "//trim(output)//"_tensCls.dat"          
     call fp%open(trim(output)//"_tensCls.dat","w")
     write(fp%unit, "(A8, 5A16)") "# ell ", "   TT  ",  "   EE  ",  " BB  ", "  TE "
     do l = lmin, lmax
        lnorm =  l*(l+1.d0)/coop_2pi*norm
        write(fp%unit, "(I8, 6E16.7)") l, tensCls(coop_index_ClTT, l)*lnorm,  tensCls(coop_index_ClEE, l)*lnorm,  tensCls(coop_index_ClBB, l)*lnorm,  tensCls(coop_index_ClTE, l)*lnorm
     enddo
     call fp%close()
     lensedCls = lensedCls + tensCls
     call fp%open(trim(output)//"_totCls.dat","w")
     write(*,*) "total Cls are saved  in "//trim(output)//"_totCls.dat"     
     write(fp%unit, "(A8, 5A16)") "# ell ", "   TT  ",  "   EE  ",  " BB  ", "  TE "
     do l=lmin, lmax
        lnorm = l*(l+1.d0)/coop_2pi*norm
        write(fp%unit, "(I5, 20E16.7)") l, lensedCls(coop_index_ClTT, l)*lnorm, lensedCls(coop_index_ClEE, l)*lnorm,  lensedCls(coop_index_ClBB, l)*lnorm,  lensedCls(coop_index_ClTE, l)*lnorm
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
