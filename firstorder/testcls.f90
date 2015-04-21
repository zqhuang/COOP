program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::fod
  type(coop_pert_object)::pert
  type(coop_file)fp
  integer,parameter::lmin=2, lmax=2600
  integer i, m, iq, ik,j, l, isrc
  COOP_REAL, dimension(coop_num_Cls, lmin:lmax)::Cls_scalar, Cls_tensor, Cls_lensed
  COOP_REAL norm
  COOP_REAL z, a, s, stau, k
  !!set cosmology
  !! call fod%Set_Planck_bestfit()
 call fod%set_standard_cosmology(Omega_b=0.0485374d0, Omega_c=0.2585497252d0, h = 0.8d0, tau_re = 0.08193d0, nu_mass_eV = 0.06d0, As = 2.2098d-9, ns = 0.968d0)
  norm = fod%cosmomc_theta()
  print*, fod%h(), 11425.8/(fod%angular_diameter_distance(1.d0/1.04d0)/fod%H0Mpc()), 3.811*fod%h()/(0.04*(1.d0+0.02*fod%HdotbyHsq(1.d0))/1.04d0)
  stop
  coop_zeta_single_slice = .true.
  !!print*, fod%zre
  !!if you want extended models
  !!
  
  norm = 2.72558**2*1.d12
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
  if(fod%index_massivenu .ne. 0)then
     call fp%open('mnu_scalCls.txt', 'w')
  else
     call fp%open('lcdm_scalCls.txt', 'w')
  endif
  do l=lmin, min(lmax, 2600)
     write(fp%unit, "(I5, 20E16.7)") l, Cls_scalar(:, l)*(l*(l+1.d0)/coop_2pi*norm), fod%clzetazeta_at_r(l, fod%distlss)*l*(l+1.)/coop_2pi*norm
  enddo
  call fp%close()
  Cls_lensed = Cls_lensed+Cls_scalar
  if(fod%index_massivenu .ne. 0)then
     call fp%open('mnu_lensedCls.txt', 'w')
  else
     call fp%open('lcdm_lensedCls.txt', 'w')
  endif
  do l=lmin, min(lmax, 2600)
     write(fp%unit, "(I5, 20E16.7)") l, Cls_lensed(1:4, l)*(l*(l+1.d0)/coop_2pi*norm)
  enddo
  call fp%close()

!!$


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
