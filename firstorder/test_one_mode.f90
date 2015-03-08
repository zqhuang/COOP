program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::fod
  COOP_INT,parameter::lmin = 2, lmax = 1000
  COOP_REAL, dimension(coop_num_Cls, lmin:lmax)::Cls_scalar, Cls_tensor, Cls_lensed
  COOP_INT l, ik
  COOP_REAL norm
  type(coop_file)fp

  norm = 2.72558**2*1.d12
  call fod%set_standard_cosmology(Omega_b=0.0485374d0, Omega_c=0.2585497252d0, h = 0.67766d0, tau_re = 0.08193d0, As = 2.2098d-9, ns = 0.968d0, nrun = 0.05d0, r = 0.d0, nt = -0.01d0, YHe = 0.248d0, Nnu = 3.d0, de_Q = 0.d0, de_tracking_n = 0.d0, de_dlnQdphi = 0.d0, de_dUdphi = -1.2d0 )
  !!***************************************************
  !! V = V0 / phi^n exp( C  * phi)
  !! n = de_tracking_n;  C = de_dUdphi
  !! V0 is determined by the condition Omega_k =0
  !!***************************************************
  !! Q = Q0 exp( A * phi)
  !! Q0 = de_Q,  A = de_dlnQdphi
  !!***************************************************

!!$!!test energy conservation
  call fod%init_source(0)  
  ik = 1
  do while(fod%source(0)%k(ik).lt. 1.d0)
     ik = ik + 1
  enddo
  call fod%compute_source_k(fod%source(0), ik, do_test_energy_conservation = .true.)
  stop
  
  
  !!compute Cl's
  call fod%compute_source(0)
  call fod%source(0)%get_All_Cls(lmin, lmax, Cls_scalar)
  call fp%open('coupled_DE_scalCls_Q0.txt', 'w')
  do l=lmin, lmax
     write(fp%unit, "(I5, 20E16.7)") l, Cls_scalar(:, l)*(l*(l+1.d0)/coop_2pi*norm), fod%clzetazeta_at_r(l, fod%distlss)*l*(l+1.)/coop_2pi*norm
  enddo
  call fp%close()
  

  
end program test
