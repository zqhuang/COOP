program simmaps
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools

  implicit none
#include "constants.h"
  COOP_INT,parameter::lmax  = 300, nsims = 100
  COOP_INT,parameter::lmin = 2, lmax_mask = 50
  COOP_INT,parameter::nside = 512
  COOP_REAL,parameter::beam_fwhm = 20.
  type(coop_healpix_maps)::map,  mask
  integer l, m, il, i, isim
  type(coop_file)::fp, fp_sim
  COOP_REAL::kernel(lmin:lmax, lmin:lmax, 4), Cl_mask(0:lmax_mask), Cl_pseudo(lmin:lmax, 6), Cl_true(lmin:lmax, 6)
  COOP_UNKNOWN_STRING, parameter::prefix = "planck15/sim"
  COOP_REAL:: Cls(lmin:lmax, 6), Cls_total(lmin:lmax, 6)
  call coop_random_init()
  Cls = 0.d0
  call fp%open_skip_comments("planck14best_lensedCls.dat")
  do l=lmin, lmax
     read(fp%unit, *) il, Cls(l,coop_TEB_index_TT), Cls(l,coop_TEB_index_EE), Cls(l,coop_TEB_index_BB), Cls(l,coop_TEB_index_TE)
     if(il.ne. l) stop "Cl file error"     
     Cls(l,:) = Cls(l,:)*(coop_2pi/l/(l+1.d0))*coop_gaussian_filter(beam_fwhm, l)**2*coop_highpass_filter(10, 20, l)**2
     if(Cls(l,coop_TEB_index_TE)**2 .gt. Cls(l,coop_TEB_index_TT)*Cls(l,coop_TEB_index_EE))then
        write(*,*) "l = ", l, " TE^2 > TT * EE"
        stop
     endif
  enddo
  call fp%close()
  call map%init(nside = nside, nmaps=3, genre="IQU", lmax = lmax)
  call mask%generate_latcut_mask(nside = nside, latitude_deg = 30.d0, depth_deg = 1.d0)
  call mask%convert2ring()
  call mask%map2alm(lmax = lmax_mask)
  cl_mask(0:lmax_mask) = mask%cl(0:lmax_mask, 1)
  call coop_pseudoCl_get_kernel_pol(lmax_mask = lmax_mask, cl_mask = cl_mask, lmin = lmin, lmax = lmax, kernel = kernel)
  Cls_total = 0.d0
  write(*,*) cls(lmax/2, 1:6)
  do isim = 1, nsims
     map%Cl = 0.
     map%Cl(lmin:lmax, 1:6) = Cls(lmin:lmax,1:6)   
     call map%simulate()
     call map%convert2ring()
     do i=1, map%nmaps
        map%map(:,i) = map%map(:,i)*mask%map(:,1)
     enddo
     call map%map2alm(lmax = lmax)
     cl_pseudo(lmin:lmax, 1:6) = map%cl(lmin:lmax, 1:6)
!     cl_true = cl_pseudo
     call coop_pseudocl2cl_pol(lmin = lmin, lmax = lmax, kernel = kernel, cl_pseudo = cl_pseudo, cl = cl_true)
     cls_total(lmin:lmax, 1:6) = cls_total(lmin:lmax, 1:6)  + cl_true
     if(mod(isim,10).eq.0)write(*,*) cls_total(lmax/2, 1:6)/isim
  enddo
  cls_total = cls_total/nsims
  call fp%open("cls_used.dat")
  call fp_sim%open("cls_sim.dat")
  do l = lmin, lmax
     write(fp%unit,"(I8, 6E16.7)") l, cls(l,:)*l*(l+1.)/coop_2pi
     write(fp_sim%unit,"(I8, 6E16.7)") l, cls_total(l, :)*l*(l+1.)/coop_2pi
  enddo
  call fp%close()
  call fp_sim%close()
end program simmaps
