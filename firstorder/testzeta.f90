program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::fod
  COOP_INT,parameter::lmax = 2500
  COOP_REAL Cls(coop_num_cls, 0:lmax)
  COOP_INT::l, isrc, i
  COOP_REAL,dimension(:,:),allocatable::cov, trans, testcov
  type(coop_file)::fp
#define DO_OUTPUT
  !!set cosmology
  call fod%Set_Planck_bestfit()
!  coop_scalar_lmax = lmax
  allocate(testcov(2,2))
  call coop_get_zeta_shells_cov(fod, 30, 2, (/0.91d0, 0.9d0 /), testcov)
  call fod%compute_source(0)
  call fod%source(0)%get_all_Cls( 0, lmax, Cls)
  call coop_set_default_zeta_r(fod%source(0), 50)

  print*, "r_max/chi_max = ", coop_zeta_r(1)/fod%source(0)%chi(1)
  allocate(cov(coop_zeta_nr, coop_zeta_nr), trans(coop_zeta_nr, fod%source(0)%nsrc))
  print*, "enter l = "
  read(*,*) l
#ifdef DO_OUTPUT
  call coop_prtsystime(.true.)
  call coop_get_zeta_trans_l(fod%source(0),  l, coop_zeta_nr, coop_zeta_r, trans)
  call coop_prtsystime()
!!$
  call coop_get_zeta_shells_cov(fod, l, coop_zeta_nr, coop_zeta_r, cov)
  trans(:,1)=  trans(:,1)*coop_zeta_dr
  trans(:,2)=  trans(:,2)*coop_zeta_dr*sqrt((l+2.d0)*(l+1.d0)*l*(l-1.d0))
  print*, dot_product(trans(:, 1), matmul(cov, trans(:, 1))), Cls(coop_index_clTT, l)
  print*, dot_product(trans(:, 2), matmul(cov, trans(:, 2))), Cls(coop_index_clEE, l)
  print*, dot_product(trans(:, 1), matmul(cov, trans(:, 2))), Cls(coop_index_clTE, l)


  call fp%open("l"//COOP_STR_OF(l)//"trans_"//COOP_STR_OF(coop_k_dense_fac)//".txt", "w")
  do i = coop_zeta_nr, 1, -1
     write(fp%unit, "(3E16.7)") coop_zeta_r(i), trans(i, :)
  enddo
  call fp%close()
#else
  do while(l.ge.2)
     call coop_get_zeta_shells_cov(fod, l, coop_zeta_nr, coop_zeta_r, cov)
     call coop_get_zeta_trans_l(fod%source(0),  l, coop_zeta_nr, coop_zeta_r, trans)
     trans(:,1)=  trans(:,1)*coop_zeta_dr
     trans(:,2)=  trans(:,2)*coop_zeta_dr*sqrt((l+2.d0)*(l+1.d0)*l*(l-1.d0))
     print*, dot_product(trans(:, 1), matmul(cov, trans(:, 1))), Cls(coop_index_clTT, l)
     print*, dot_product(trans(:, 2), matmul(cov, trans(:, 2))), Cls(coop_index_clEE, l)
     print*, dot_product(trans(:, 1), matmul(cov, trans(:, 2))), Cls(coop_index_clTE, l)
     print*, "enter l = "
     read(*,*) l
  enddo
#endif
!!$  call coop_generate_3Dzeta(fod, "zeta3D.dat", lmax, fod%source(0)%ntau, fod%source(0)%chi)
!!$
!!$  call coop_load_2Dzeta("zeta3D.dat", lmax, fod%source(0)%ntau, 30, alms)

end program test
