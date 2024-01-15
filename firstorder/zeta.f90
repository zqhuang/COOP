program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::fod
  COOP_INT::l, isrc, i
  COOP_REAL,dimension(:,:),allocatable::cov, trans, testcov
  type(coop_file)::fp
  !!set cosmology
  call fod%Set_Planck_bestfit()
  allocate(testcov(2,2))
  call coop_get_zeta_shells_cov(fod, 30, 2, (/0.91d0, 0.9d0 /), testcov)
  call fod%compute_source(0)
  call coop_set_default_zeta_r(fod%source(0), 50)

  print*, "r_max/chi_max = ", coop_zeta_r(1)/fod%source(0)%chi(1)
  allocate(trans(coop_zeta_nr, fod%source(0)%nsrc))
  do l = 2, 100
     write(*,*) "doing l = ", l
     call coop_get_zeta_trans_l(fod%source(0),  l, coop_zeta_nr, coop_zeta_r, trans)
     trans(:,1)=  trans(:,1)*coop_zeta_dr
     call fp%open("trans_"//COOP_STR_OF(l)//".dat","u")
     write(fp%unit) trans
     call fp%close()
  enddo
!!$  call coop_generate_3Dzeta(fod, "zeta3D.dat", lmax, fod%source(0)%ntau, fod%source(0)%chi)
!!$
!!$  call coop_load_2Dzeta("zeta3D.dat", lmax, fod%source(0)%ntau, 30, alms)

end program test
