program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::fod
  type(coop_pert_object)::pert
  type(coop_file)fp
  integer,parameter::lmin=2, lmax=2700
  integer i, m, iq, ik,j, l, isrc
  COOP_REAL, dimension(coop_num_Cls, lmin:lmax)::Cls_scalar, Cls_tensor, Cls_lensed
  COOP_REAL norm, hub, lnorm
  COOP_REAL z, a, s, stau, k
  !!set cosmology
  !! call fod%Set_Planck_bestfit()
  hub = 0.676d0
 call fod%set_standard_cosmology(Omega_b = 0.022d0/hub**2, omega_c = 0.12d0/hub**2, h = hub, tau_re = 0.08d0, As = 2.21979d-9, ns = 0.96d0)
!  coop_zeta_single_slice = .true.
  !!print*, fod%zre
  !!if you want extended models
  !!
  
  norm = fod%Tcmb()**2*1.d12
  !!compute the scalar Cl's
  call coop_prtSystime(.true.)
  call fod%compute_source(0)
  call coop_prtSystime()
!!$

  call coop_prtSystime(.true.)
  call fod%source(0)%get_All_Cls(lmin, lmax, Cls_scalar)
  call coop_prtSystime()
  call coop_prtSystime(.true.)
  call coop_get_lensing_Cls(lmin, lmax, Cls_Scalar, Cls_lensed)
  call coop_prtSystime()
  call fp%open('std6_cls.txt', 'w')
  do l=lmin, min(lmax, 2500)
     lnorm = l*(l+1.d0)/coop_2pi*norm
     write(fp%unit, "(I5, 20E16.7)") l, Cls_scalar(coop_index_ClTT, l)*lnorm, Cls_scalar(coop_index_ClEE, l)*lnorm, Cls_scalar(coop_index_ClTE, l)*lnorm, Cls_scalar(coop_index_ClLenLen, l)*norm*dble(l)**4, Cls_scalar(coop_index_ClTLen, l)*norm*dble(l)**4
  enddo
  call fp%close()
!  Cls_lensed = Cls_lensed+Cls_scalar

!!$
  stop

  !!compute the tensor Cl's
  call coop_prtSystime(.true.)
  call fod%compute_source(2)
  call coop_prtSystime()

  call coop_prtSystime(.true.)
  call fod%source(2)%get_All_Cls(lmin,lmax, Cls_tensor)
  call coop_prtSystime()
  if(fod%index_massivenu .ne. 0)then
     call fp%open('mnu_tensCls.txt', 'w')
  else
     call fp%open('lcdm_tensCls.txt', 'w')
  endif
  do l=lmin, min(lmax, 1500)
     write(fp%unit, "(I5, 20E16.7)") l, Cls_tensor(1:4, l)*(l*(l+1.d0)/coop_2pi*norm)
  enddo
  call fp%close()

  Cls_lensed = Cls_lensed + Cls_tensor

  if(fod%index_massivenu .ne. 0)then
     call fp%open('mnu_totCls.txt', 'w')
  else
     call fp%open('lcdm_totCls.txt', 'w')
  endif
  do l=lmin, min(lmax, 2600)
     write(fp%unit, "(I5, 20E16.7)") l,  Cls_lensed(1:4, l)*(l*(l+1.d0)/coop_2pi*norm)
  enddo
  call fp%close()


end program test
