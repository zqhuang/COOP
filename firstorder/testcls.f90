program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::cosmology
  type(coop_file)fp
  integer,parameter::lmin=2, lmax=2800, nk=256
  integer i, l, ell_now
  COOP_REAL, dimension(coop_num_Cls, lmin:lmax)::Cls_scalar, Cls_tensor, Cls_lensed
  COOP_REAL norm, hub, lnorm
  COOP_REAL k(nk), matterPk(nk), khMpc(nk)
  COOP_REAL,dimension(:),allocatable::intpp 
  !!set cosmology
  !! call cosmology%Set_Planck_bestfit()
  hub = 0.676d0
 call cosmology%set_standard_cosmology(Omega_b = 0.022d0/hub**2, omega_c = 0.12d0/hub**2, h = hub, tau_re = 0.08d0, As = 2.21979d-9, ns = 0.96d0)
!  coop_zeta_single_slice = .true.
  !!print*, cosmology%zre
  !!if you want extended models
  !!
  norm = cosmology%Tcmb()**2*1.d12
  !!compute the scalar Cl's
  call coop_prtSystime(.true.)
  call cosmology%compute_source(0)
  call coop_prtSystime()

  call cosmology%Set_Planck_bestfit()
  call coop_prtSystime(.true.)
  call cosmology%compute_source(0)
  call coop_prtSystime()

  
  call coop_prtSystime(.true.)
  call cosmology%source(0)%get_All_Cls(lmin, lmax, Cls_scalar)
  call coop_prtSystime()
  call coop_prtSystime(.true.)
  call coop_get_lensing_Cls(lmin, lmax, Cls_Scalar, Cls_lensed)
  call coop_prtSystime()
  call fp%open('std6_scalCls.txt', 'w')
  do l=lmin, min(lmax, 2500)
     lnorm = l*(l+1.d0)/coop_2pi*norm
     write(fp%unit, "(I5, 20E16.7)") l, Cls_scalar(coop_index_ClTT, l)*lnorm, Cls_scalar(coop_index_ClEE, l)*lnorm, Cls_scalar(coop_index_ClTE, l)*lnorm, Cls_scalar(coop_index_ClLenLen, l)*norm*dble(l)**4, Cls_scalar(coop_index_ClTLen, l)*norm*dble(l)**3
  enddo
  call fp%close()
  Cls_lensed = Cls_lensed+Cls_scalar
  call fp%open('std6_lensedCls.txt', 'w')
  do l=lmin, min(lmax, 2500)
     lnorm = l*(l+1.d0)/coop_2pi*norm
     write(fp%unit, "(I5, 20E16.7)") l, Cls_lensed(coop_index_ClTT, l)*lnorm, Cls_lensed(coop_index_ClEE, l)*lnorm,  Cls_lensed(coop_index_ClBB, l)*lnorm,  Cls_lensed(coop_index_ClTE, l)*lnorm
  enddo
  call fp%close()

  call coop_set_uniform(nk, k, 10.d0, cosmology%source(0)%kmax, logscale = .true.)
  khMpc = k * cosmology%H0Mpc()/cosmology%h()  !!k/H0 * (H0 * Mpc) / h = k in unit of h Mpc^{-1}
  call cosmology%get_matter_power(z=0.d0, nk = nk, k = k, Pk = matterPk)  !!this returns k^3 |\delta_k|^2 /(2pi^2)
  matterPk = matterPk * (2.d0*coop_pi**2)/khMpc**3   
  call fp%open('std6_matterpower.txt', 'w')
  do i = 1, nk
     write(fp%unit, "(2E16.7)") khMpc(i),  matterPk(i)
  enddo
  call fp%close()
  
!!$
  stop

  !!compute the tensor Cl's
  call coop_prtSystime(.true.)
  call cosmology%compute_source(2)
  call coop_prtSystime()

  call coop_prtSystime(.true.)
  call cosmology%source(2)%get_All_Cls(lmin,lmax, Cls_tensor)
  call coop_prtSystime()
  call fp%open('std6_tensCls.txt', 'w')
  do l=lmin, min(lmax, 1500)
     write(fp%unit, "(I5, 20E16.7)") l, Cls_tensor(1:4, l)*(l*(l+1.d0)/coop_2pi*norm)
  enddo
  call fp%close()
  
  Cls_lensed = Cls_lensed + Cls_tensor

  call fp%open('std6_totCls.txt', 'w')
  do l=lmin, min(lmax, 2600)
     write(fp%unit, "(I5, 20E16.7)") l,  Cls_lensed(1:4, l)*(l*(l+1.d0)/coop_2pi*norm)
  enddo
  call fp%close()


contains


  function intphiphi(chi)
    COOP_REAL chi, intphiphi, phi, k
    if(chi .lt. 1.d-4)then
       intphiphi = 0.d0
       return
    endif
    k = ell_now/chi
    phi = cosmology%source(0)%interpolate_one(k, chi, coop_index_source_Len)
    intphiphi = cosmology%psofk(ell_now/chi)*chi*phi**2
  end function intphiphi


end program test
