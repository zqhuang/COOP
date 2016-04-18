program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  use coop_fitsio_mod
  implicit none
#include "constants.h"
  COOP_INT,parameter::lmax = 500, lmin = 40, lmax_mask = 200, delta_l = 15
  type(coop_healpix_maps)::hm1, hm2, mask, hm
  COOP_REAL,parameter::nuwidth = 0.1d0, ewidth = 0.03d0
  COOP_REAL::nucut, ecut
  COOP_INT::i, l, lp, lfilter
  COOP_REAL,dimension(:,:,:),allocatable::kernel
  COOP_REAL,dimension(:),allocatable::Cl_mask
  COOP_REAL,dimension(:,:),allocatable::Cl1_pseudo, Cl2_pseudo, ClCross_pseudo, Cl1, Cl2, ClCross
  COOP_REAL,dimension(:),allocatable::ells, binned_Cls, ClEE, ClBB
  COOP_REAL::rms1, rms2, mean1, mean2, sm, weight, mean, rms
  type(coop_asy)::fig
  COOP_STRING::prefix
  call coop_get_command_line_argument(key = "l", arg = lfilter)
!!$  call coop_get_command_line_argument(key = "nu", arg = nucut)
!!$  call coop_get_command_line_argument(key = "e", arg = ecut)
!!$  call hm%read("dust/dust_i_30a_hp_20_40_n1024_conv_TQUL.fits")
  call hm1%read("dust/dust_iqu_10a_hp_20_40_n1024_hm1.fits")
  call hm2%read("dust/dust_iqu_10a_hp_20_40_n1024_hm2.fits")
!!$  call fig%open("EB_nucut"//COOP_NICESTR_OF(nucut)//"_ecut"//COOP_NICESTR_OF(ecut)//".txt")
  call fig%open("EB_l"//COOP_STR_OF(lfilter)//".txt")
  prefix = ""
  call fig%init(xlabel = "$\ell$", ylabel = "$\ell(\ell+1)C_\ell/(2\pi)$")
  call fig%line(xstart = dble(lmin), ystart = 0.d0, xend = dble(lmax), yend = 0.d0, color = "black", linewidth = 1.5)

  call mask%read("dust/mask_l"//COOP_STR_OF(lfilter)//".fits")
  allocate(kernel(lmin:lmax, lmin:lmax, 4), Cl_mask(0:lmax_mask), CL1_pseudo(lmin:lmax, 6), Cl1(lmin:lmax, 6),  CL2_pseudo(lmin:lmax, 6), Cl2(lmin:lmax, 6),  CLCross_pseudo(lmin:lmax, 6), ClCross(lmin:lmax, 6), ells(lmin:lmax), binned_Cls(lmin:lmax), ClEE(lmin:lmax), ClBB(lmin:lmax))

!!$  hm%map(:,1) = log(max(hm%map(:,1), 1.d-5))
!!$  call hm%apply_mask(mask)
  call hm1%apply_mask(mask)
  call hm2%apply_mask(mask)

!!$  !!do threshold cut, produce new mask 
!!$  sm = sum(mask%map(:,1))
!!$  mean = sum(hm%map(:,1))/sm
!!$  rms = sqrt(sum((hm%map(:,1)-mean)**2*mask%map(:,1))/sm)
!!$  print*, "mask fsky = "//COOP_STR_OF(sum(mask%map(:,1))/mask%npix)
!!$  if(ecut .gt. 0.d0)then
!!$     mask%map(:,1) = mask%map(:,1)*(1.d0 - (1.d0+tanh(((hm%map(:,1)-mean)/rms-nucut)/nuwidth))/2.d0 *(1.d0+tanh( (sqrt(hm%map(:,2)**2+hm%map(:,3)**2)/max(abs(hm%map(:,4)), sqrt(hm%map(:,2)**2+hm%map(:,3)**2)+1.d-99) - ecut)/ewidth) )/2.d0)
!!$  else
!!$     mask%map(:,1) = mask%map(:,1)*(1.d0 - (1.d0+tanh(((hm%map(:,1)-mean)/rms-nucut)/nuwidth))/2.d0)
!!$  endif
!!$
!!$  print*, "after threshold cut, mask fsky = "//COOP_STR_OF(sum(mask%map(:,1))/mask%npix)
!!$  call mask%write("mask.fits")
!!$  call hm1%apply_mask(mask)
!!$  call hm2%apply_mask(mask)
  !!do threshold cut
!!$  sm = sum(mask%map(:,1))
!!$  mean1  = sum(hm1%map(:,1))/sm
!!$  rms1 = sqrt(sum((hm1%map(:,1)-mean1)**2*mask%map(:,1))/sm)
!!$  mean2  = sum(hm2%map(:,1))/sm
!!$  rms2 = sqrt(sum((hm2%map(:,1)-mean2)**2*mask%map(:,1))/sm)
!!$  hm1%map(:,1) = hm1%map(:,1)*(1.+tanh( ((hm1%map(:,1)-mean1)/rms1 - nucut )/nuwidth))/2.*mask%map(:,1)
!!$  hm2%map(:,1) = hm2%map(:,1)*(1.+tanh( ((hm2%map(:,1)-mean2)/rms2 - nucut )/nuwidth))/2.*mask%map(:,1)
!!$  call hm1%apply_mask(mask)
!!$  call hm2%apply_mask(mask)

  call hm1%map2alm(lmax = lmax)
  call hm2%map2alm(lmax = lmax)
  call hm1%get_cls()
  Cl1_pseudo(lmin:lmax, 1:6) = hm1%cl(lmin:lmax, 1:6)
  call hm2%get_cls()
  Cl2_pseudo(lmin:lmax, 1:6) = hm2%cl(lmin:lmax, 1:6)
  call hm1%get_cls(cross = hm2)
  ClCross_pseudo(lmin:lmax, 1:6) = hm1%cl(lmin:lmax, 1:6)

  call mask%map2alm(lmax = lmax_mask)
  Cl_mask(0:lmax_mask) = mask%cl(0:lmax_mask,1)
  call coop_pseudoCl_get_kernel_pol(lmax_mask = lmax_mask, Cl_mask = Cl_mask, lmin = lmin, lmax = lmax, kernel = kernel)

  call coop_pseudoCl2Cl_pol(lmin = lmin, lmax = lmax, Cl_pseudo = Cl1_pseudo, kernel = kernel, Cl = Cl1, smooth = .false.)
  call coop_pseudoCl2Cl_pol(lmin = lmin, lmax = lmax, Cl_pseudo = Cl2_pseudo, kernel = kernel, Cl = Cl2, smooth = .false.)
  call coop_pseudoCl2Cl_pol(lmin = lmin, lmax = lmax, Cl_pseudo = ClCross_pseudo, kernel = kernel, Cl = ClCross, smooth = .false.)



  ClEE = ClCross(:, coop_TEB_index_EE) - Cl1(:, coop_TEB_index_TE)*Cl2(:, coop_TEB_index_TE)*ClCross(:, coop_TEB_index_TT)/Cl1(:, coop_TEB_index_TT)/Cl2(:, coop_TEB_index_TT)
  ClBB = ClCross(:, coop_TEB_index_BB) - Cl1(:, coop_TEB_index_TB)*Cl2(:, coop_TEB_index_TB)*ClCross(:, coop_TEB_index_TT)/Cl1(:, coop_TEB_index_TT)/Cl2(:, coop_TEB_index_TT)

  print*, "unconstrained EE"
  do l = lmin, lmax
     ells(l) = l
     binned_Cls(l) = 0.d0
     weight = 0.d0
     do lp = max(lmin, l-3*delta_l), min(l+3*delta_l, lmax)
        binned_Cls(l) =  binned_Cls(l) + lp*(lp+1.d0)*ClEE(lp)/coop_2pi * exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
        weight = weight + exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
     enddo
     binned_Cls(l) = binned_Cls(l)/weight
     if(mod(l, 50).eq.0) print*, l, binned_Cls(l)
  enddo
  call fig%plot(ells(lmin+delta_l:lmax-delta_l), binned_Cls(lmin+delta_l:lmax-delta_l), color = "RGB:200:50:50", legend="unconstrained EE", linewidth = 2., linetype = "solid")
  print*, "unconstrained BB"
  do l = lmin, lmax
     ells(l) = l
     binned_Cls(l) = 0.d0
     weight = 0.d0
     do lp = max(lmin, l-3*delta_l), min(l+3*delta_l, lmax)
        binned_Cls(l) =  binned_Cls(l) + lp*(lp+1.d0)*ClBB(lp)/coop_2pi * exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
        weight = weight + exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
     enddo
     binned_Cls(l) = binned_Cls(l)/weight
     if(mod(l, 50).eq.0) print*, l, binned_Cls(l)
  enddo
  call fig%plot(ells(lmin+delta_l:lmax-delta_l), binned_Cls(lmin+delta_l:lmax-delta_l), color = "RGB:50:50:200", linewidth = 2.5, linetype = "dotted", legend="unconstrained BB")


  print*, "EE"
  do l = lmin, lmax
     ells(l) = l
     binned_Cls(l) = 0.d0
     weight = 0.d0
     do lp = max(lmin, l-3*delta_l), min(l+3*delta_l, lmax)
        binned_Cls(l) =  binned_Cls(l) + lp*(lp+1.d0)*ClCross(lp,coop_TEB_index_EE)/coop_2pi * exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
        weight = weight + exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
     enddo
     binned_Cls(l) = binned_Cls(l)/weight
     if(mod(l, 50).eq.0) print*, l, binned_Cls(l)
  enddo
  call fig%plot(ells(lmin+delta_l:lmax-delta_l), binned_Cls(lmin+delta_l:lmax-delta_l), color = "RGB:50:200:50", legend="EE", linewidth = 1., linetype = "solid")
  print*, "BB"
  do l = lmin, lmax
     ells(l) = l
     binned_Cls(l) = 0.d0
     weight = 0.d0
     do lp = max(lmin, l-3*delta_l), min(l+3*delta_l, lmax)
        binned_Cls(l) =  binned_Cls(l) + lp*(lp+1.d0)*ClCross(lp, coop_TEB_index_BB)/coop_2pi * exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
        weight = weight + exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
     enddo
     binned_Cls(l) = binned_Cls(l)/weight
     if(mod(l, 50).eq.0) print*, l, binned_Cls(l)
  enddo
  call fig%plot(ells(lmin+delta_l:lmax-delta_l), binned_Cls(lmin+delta_l:lmax-delta_l), color = "cyan", linewidth = 1.5, linetype = "dotted", legend="BB")

!!$  write(*,*) "EB"
!!$  do l = lmin, lmax
!!$     ells(l) = l
!!$     binned_Cls(l) = 0.d0
!!$     weight = 0.d0
!!$     do lp = max(lmin, l-3*delta_l), min(l+3*delta_l, lmax)
!!$        binned_Cls(l) =  binned_Cls(l) + lp*(lp+1.d0)*ClCross(lp, coop_TEB_index_EB)/coop_2pi * exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
!!$        weight = weight + exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
!!$     enddo
!!$     binned_Cls(l) = binned_Cls(l)/weight
!!$     if(mod(l, 100).eq.0) print*, l, binned_Cls(l), Cl_pseudo(l, coop_TEB_index_EB)*l*(l+1.)/coop_2pi
!!$  enddo
!!$  call fig%plot(ells(lmin+delta_l:lmax-delta_l), binned_Cls(lmin+delta_l:lmax-delta_l), color = "RGB:50:200:50", legend=trim(prefix)//"EB", linewidth = 1.5, linetype = "solid")
  call fig%legend(0.03, 0.26)
  call fig%close()
  
!!$     
!!$
!!$
!!$  write(*,*) "TT"
!!$  do l = lmin, lmax
!!$     ells(l) = l
!!$     binned_Cls(l) = 0.d0
!!$     weight = 0.d0
!!$     do lp = max(lmin, l-3*delta_l), min(l+3*delta_l, lmax)
!!$        binned_Cls(l) =  binned_Cls(l) + lp*(lp+1.d0)*Cl(lp, coop_TEB_index_TT)/coop_2pi * exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
!!$        weight = weight + exp(-dble(l-lp)**2/(dble(delta_l)**2*2.))
!!$     enddo
!!$     binned_Cls(l) = binned_Cls(l)/weight
!!$     if(mod(l, 100).eq.0) print*, l, binned_Cls(l), Cl_pseudo(l, coop_TEB_index_TT)*l*(l+1.)/coop_2pi
!!$  enddo

!!# Jysr_to_muKRJ = c**2 / 2. / K_b / (freq * 1e9) **2 * 1e-26 * 1e6 #  [uK_RJ/(Jy/sr)]


end program test
