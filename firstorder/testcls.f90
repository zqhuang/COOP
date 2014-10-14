program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::fod
  type(coop_pert_object)::pert
  type(coop_file)fp
  integer,parameter::lmin=2, lmax=2600
  integer i, m, iq, ik,j, l, isrc
  COOP_REAL, dimension(coop_num_Cls, lmin:lmax)::Cls_scalar, Cls_tensor
  COOP_REAL norm
  COOP_REAL z, a, s, stau, k

  !!set cosmology
  call fod%Set_Planck_bestfit_with_r(r=0.2d0, Omega_nu = 0.020176d0)
  !!print*, fod%zre
  !!if you want extended models
  !!call fod%set_standard_cosmology(Omega_b=0.047d0, Omega_c=0.26d0, h = 0.68d0, tau_re = 0.08d0, nu_mass_eV = 0.06d0, As = 2.15d-9, ns = 0.962d0, nrun = -0.01d0, r = 0.2d0, nt = -0.01d0, YHe = 0.25d0, Nnu = 3.d0)
  norm = 2.72558**2*1.d12

  !!compute the scalar Cl's
  call coop_prtSystime(.true.)
  call fod%compute_source(0)
  call coop_prtSystime()
!!$
  call coop_prtSystime(.true.)
  call fod%source(0)%get_All_Cls(lmin, lmax, Cls_scalar)
  call coop_prtSystime()
  if(fod%index_massivenu .ne. 0)then
     call fp%open('mnu_scalCls.txt', 'w')
  else
     call fp%open('lcdm_scalCls.txt', 'w')
  endif
  do l=lmin, min(lmax, 2600)
     write(fp%unit, "(I5, 20E16.7)") l, Cls_scalar(1, l)*(l*(l+1.d0)/coop_2pi*norm), Cls_scalar(2, l)*(l*(l+1.d0)/coop_2pi*norm), Cls_scalar(4, l)*(l*(l+1.d0)/coop_2pi*norm), Cls_scalar(5, l)*real(l)**4*norm, Cls_scalar(6, l)*real(l)**3*norm
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


end program test
