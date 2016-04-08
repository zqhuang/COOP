program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  use coop_fitsio_mod
  implicit none
#include "constants.h"
  logical::do_log = .false.
  COOP_INT,parameter::lmax = 500, lmin = 30, lmax_mask = 200, delta_l = 8
  type(coop_healpix_maps)::hm1, hm2, mask
  COOP_REAL, parameter::nucut = 1.d0, nuwidth = 0.1d0
  COOP_INT::i, l, lp
  COOP_REAL,dimension(:,:,:),allocatable::kernel
  COOP_REAL,dimension(:),allocatable::Cl_mask
  COOP_REAL,dimension(:,:),allocatable::Cl_pseudo, Cl
  COOP_REAL,dimension(:),allocatable::ells, binned_Cls
  COOP_REAL::rms, weight
  type(coop_asy)::fig
  COOP_STRING::prefix
  call coop_get_command_line_argument(key = "log", arg = do_log, default = .false.)
  if(do_log)then
     call hm1%read("planck15/dust_iqu_10a_n1024_hm1_ln.fits")
     call hm2%read("planck15/dust_iqu_10a_n1024_hm2_ln.fits")
     call fig%open("lognorm_EB.txt")
     prefix = "log-norm "
  else
     call hm1%read("planck15/dust_iqu_10a_n1024_hm1.fits")
     call hm2%read("planck15/dust_iqu_10a_n1024_hm2.fits")
     call fig%open("EB.txt")
     prefix = ""
  endif
  call fig%init(xlabel = "$\ell$", ylabel = "$\ell(\ell+1)C_\ell/(2\pi)$", ymin = 0.)
  call fig%line(xstart = dble(lmin), ystart = 0.d0, xend = dble(lmax), yend = 0.d0, color = "black", linewidth = 1.5)

  call mask%read("planck15/mask_lat30_n1024.fits")
  allocate(kernel(lmin:lmax, lmin:lmax, 4), Cl_mask(0:lmax_mask), CL_pseudo(lmin:lmax, 6), Cl(lmin:lmax, 6), ells(lmin:lmax), binned_Cls(lmin:lmax))
  call mask%map2alm(lmax = lmax_mask)
  Cl_mask(0:lmax_mask) = mask%cl(0:lmax_mask,1)
  call coop_pseudoCl_get_kernel_pol(lmax_mask = lmax_mask, Cl_mask = Cl_mask, lmin = lmin, lmax = lmax, kernel = kernel)

  call hm1%apply_mask(mask)
  call hm1%map2alm(lmax = lmax)
  call hm2%map2alm(lmax = lmax)
  call hm1%get_cls(cross = hm2)
  Cl_pseudo(lmin:lmax, 1:6) = hm1%cl(lmin:lmax, 1:6)
  call coop_pseudoCl2Cl_pol(lmin = lmin, lmax = lmax, Cl_pseudo = Cl_pseudo, kernel = kernel, Cl = Cl, smooth = .false.)
  print*, "EE"
  do l = lmin, lmax
     ells(l) = l
     binned_Cls(l) = 0.d0
     weight = 0.d0
     do lp = max(lmin, l-3*delta_l), min(l+3*delta_l, lmax)
        binned_Cls(l) =  binned_Cls(l) + lp*(lp+1.d0)*Cl(lp, coop_TEB_index_EE)/coop_2pi * exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
        weight = weight + exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
     enddo
     binned_Cls(l) = binned_Cls(l)/weight
     if(mod(l, 100).eq.0) print*, l, binned_Cls(l), Cl_pseudo(l, coop_TEB_index_EE)*l*(l+1.)/coop_2pi
  enddo
  call fig%plot(ells(lmin+delta_l:lmax-delta_l), binned_Cls(lmin+delta_l:lmax-delta_l), color = "RGB:200:50:50", legend=trim(prefix)//"EE", linewidth = 2., linetype = "solid")
  print*, "BB"
  do l = lmin, lmax
     ells(l) = l
     binned_Cls(l) = 0.d0
     weight = 0.d0
     do lp = max(lmin, l-3*delta_l), min(l+3*delta_l, lmax)
        binned_Cls(l) =  binned_Cls(l) + lp*(lp+1.d0)*Cl(lp, coop_TEB_index_BB)/coop_2pi * exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
        weight = weight + exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
     enddo
     binned_Cls(l) = binned_Cls(l)/weight
     if(mod(l, 100).eq.0) print*, l, binned_Cls(l), Cl_pseudo(l, coop_TEB_index_BB)*l*(l+1.)/coop_2pi
  enddo
  call fig%plot(ells(lmin+delta_l:lmax-delta_l), binned_Cls(lmin+delta_l:lmax-delta_l), color = "RGB:50:50:200", linewidth = 2.5, linetype = "dotted", legend=trim(prefix)//"BB")
  write(*,*) "EB"
  do l = lmin, lmax
     ells(l) = l
     binned_Cls(l) = 0.d0
     weight = 0.d0
     do lp = max(lmin, l-3*delta_l), min(l+3*delta_l, lmax)
        binned_Cls(l) =  binned_Cls(l) + lp*(lp+1.d0)*Cl(lp, coop_TEB_index_EB)/coop_2pi * exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
        weight = weight + exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
     enddo
     binned_Cls(l) = binned_Cls(l)/weight
     if(mod(l, 100).eq.0) print*, l, binned_Cls(l), Cl_pseudo(l, coop_TEB_index_EB)*l*(l+1.)/coop_2pi
  enddo
  call fig%plot(ells(lmin+delta_l:lmax-delta_l), binned_Cls(lmin+delta_l:lmax-delta_l), color = "RGB:50:200:50", legend=trim(prefix)//"EB", linewidth = 1.5, linetype = "solid")
  call fig%legend(0.1, 0.9)
  call fig%close()

     


  write(*,*) "TT"
  do l = lmin, lmax
     ells(l) = l
     binned_Cls(l) = 0.d0
     weight = 0.d0
     do lp = max(lmin, l-3*delta_l), min(l+3*delta_l, lmax)
        binned_Cls(l) =  binned_Cls(l) + lp*(lp+1.d0)*Cl(lp, coop_TEB_index_TT)/coop_2pi * exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
        weight = weight + exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
     enddo
     binned_Cls(l) = binned_Cls(l)/weight
     if(mod(l, 100).eq.0) print*, l, binned_Cls(l), Cl_pseudo(l, coop_TEB_index_TT)*l*(l+1.)/coop_2pi
  enddo

!!# Jysr_to_muKRJ = c**2 / 2. / K_b / (freq * 1e9) **2 * 1e-26 * 1e6 #  [uK_RJ/(Jy/sr)]


end program test
