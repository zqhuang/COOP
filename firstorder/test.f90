program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::fod
  type(coop_pert_object)::pert
  type(coop_file)fp
  integer,parameter::lmin=2, lmax=1000
  integer i, m, iq, ik,j, l, isrc
  COOP_REAL, dimension(coop_num_Cls, lmin:lmax)::Cls_scalar, Cls_tensor
  COOP_REAL norm
  COOP_REAL z, a, s, stau, k

  !!set cosmology
  call fod%Set_Planck_bestfit_with_r(r=0.2d0)
    !if you want extended models
!  call fod%set_standard_cosmology(Omega_b=0.047d0, Omega_c=0.26d0, h = 0.68d0, tau_re = 0.08d0, nu_mass_eV = 0.06d0, As = 2.15d-9, ns = 0.962d0, nrun = -0.01d0, r = 0.2d0, nt = -0.01d0, YHe = 0.25d0, Nnu = 3.d0)
  norm = 2.726**2*1.d12

!!$  m = 0
!!$  !!compute the scalar source
!!$  call coop_prtSystime(.true.)
!!$  call fod%compute_source(m)
!!$  call coop_prtSystime()
!!$  ik = 1
!!$  call fp%open("source_p1.txt")
!!$  do i= 1, fod%source(m)%ntau
!!$     write(fp%unit,"(14E16.7)") fod%source(m)%tau(i), fod%source(m)%saux(1, ik, i), fod%source(m)%s(1, ik, i)
!!$  enddo
!!$  call fp%close()
!!$  stop
!!$
!!$  call fp%open("source10.txt")
!!$  k = 4.80146
!!$  do i= 1, fod%source(m)%ntau
!!$     write(fp%unit,"(14E16.7)") fod%source(m)%tau(i), fod%source(m)%interpolate(k, fod%source(m)%chi(i))
!!$  enddo
!!$  call fp%close()
!!$
!!$
!!$  call fp%open("source30.txt")
!!$  k = 140.23
!!$  do i= 1, fod%source(m)%ntau
!!$     write(fp%unit,"(14E16.7)") fod%source(m)%tau(i), fod%source(m)%interpolate(k, fod%source(m)%chi(i))
!!$  enddo
!!$  call fp%close()

!!$
!!$  ik = 1
!!$  do while(fod%source(2)%k(ik).lt. 4.8 .and. ik.lt.fod%source(2)%nk)
!!$     ik = ik+1
!!$  enddo
!!$  print*, ik, fod%source(2)%k(ik)
!!$     
!!$  call fp%open("tensor_source_"//COOP_STR_OF(ik)//".txt", "w")
!!$  do i=1, fod%source(2)%ntau
!!$     write(fp%unit, "(20E16.7)") fod%source(2)%tau(i), fod%source(2)%s(:, ik, i)
!!$  enddo
!!$  call fp%close()
!!$
!!$
!!$  stop

  !!compute the scalar Cl's
  call coop_prtSystime(.true.)
  call fod%compute_source(0)
  call coop_prtSystime()
!!$
  call coop_prtSystime(.true.)
  call fod%source(0)%get_All_Cls(lmin, lmax, Cls_scalar)
  call coop_prtSystime()

  call fp%open('Cls_scalar.txt', 'w')
  do l=lmin, lmax
     write(fp%unit, "(I5, 20E16.7)") l, Cls_scalar(:, l)*(l*(l+1.d0)/coop_2pi*norm)
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

  call fp%open('Cls_tensor.txt', 'w')
  do l=lmin, lmax
     write(fp%unit, "(I5, 20E16.7)") l, Cls_tensor(:, l)*(l*(l+1.d0)/coop_2pi*norm)
  enddo
  call fp%close()


contains

  function trans_int(chi) result(t)
    COOP_REAL chi, t, s
    t = coop_jl(l, k*chi)*fod%source(m)%interpolate_one(k, chi, isrc)
  end function trans_int
  

end program test
