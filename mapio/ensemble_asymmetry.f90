program test
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
  COOP_UNKNOWN_STRING, parameter::prefix = "predx11"
  COOP_UNKNOWN_STRING, parameter::color_table = "Planck"
  COOP_UNKNOWN_STRING, parameter::spot_type = "Tmax"
  COOP_UNKNOWN_STRING, parameter::stack_type = "T"
  COOP_REAL, parameter::threshold = 0
  COOP_INT, parameter::mmax = 0
  COOP_INT, parameter::n = 30
  COOP_REAL, parameter::dr = 10.*coop_SI_arcmin
  COOP_INT, parameter:: imap = 1
  integer,parameter::n_sim = 1000
  COOP_STRING::fmt, fmtscreen
  COOP_UNKNOWN_STRING, parameter::resol = "512"
  COOP_UNKNOWN_STRING, parameter::map_file = prefix//"/"//prefix//"_iqu"//resol//".fits"
  COOP_UNKNOWN_STRING, parameter::imask_file  = prefix//"/"//prefix//"_imask"//resol//".fits"
  COOP_UNKNOWN_STRING, parameter::polmask_file  = prefix//"/"//prefix//"_polmask"//resol//".fits"
  COOP_UNKNOWN_STRING, parameter::fr_file = prefix//"_"//stack_type//"_sim_Fr"//resol//".dat"
  type(coop_healpix_patch)::patch_s, patch_n
  type(coop_healpix_maps)::map, sim, tmp
  type(coop_healpix_maps)::imask, polmask
  type(coop_healpix_maps)::stack_mask, spots_mask
  COOP_REAL frmean(0:n), cov(0:n,0: n), diff(0:n), chisq, err(0:n),  hdir(2)
  COOP_REAL psi(0:n, 0:n), eig(0:n)
  COOP_INT :: nspots, nlines, il, i, j, step, weight
  type(coop_asy)::fig
  type(coop_file)::fp
  type(coop_list_integer)::listpix
  type(coop_list_real)::listangle
  type(coop_file) fpsim
  integer,parameter::ncols_used = 5
  COOP_REAL, parameter::epsilon = 1.d-8
  COOP_INT nmaps_wanted, ismall, ilarge

  !!read mask and map
  call imask%read(imask_file, nmaps_wanted = 1)
  call polmask%read(polmask_file, nmaps_wanted = 1)
  nmaps_wanted = 1
  select case(trim(stack_type))
  case("I","T")
     stack_mask = imask
  case default
     stack_mask = polmask
     nmaps_wanted = 3
  end select
  select case(trim(spot_type))
  case("Tmax", "Tmin", "Tmax_QTUTOrient", "Tmin_QTUTOrient")
     spots_mask = imask
  case default
     spots_mask = polmask
     nmaps_wanted = 3
  end select

  call map%read(map_file, nmaps_wanted = nmaps_wanted)

  if(nmaps_wanted .eq. 3)then
     call map%get_fullCls(imask, polmask)
  else
     call map%get_fullCls(imask)
  endif
  call coop_healpix_lb2ang(l_deg = 226.d0, b_deg = -17.d0, theta = hdir(1), phi = hdir(2))

  sim = map
  tmp = sim

  call coop_random_init()

  frmean = 0.
  cov = 0.
  fmt = "("//trim(coop_num2str(n+1))//"G16.7)"
  fmtscreen = "("//trim(coop_num2str(n+1))//"F7.2)"
  weight = 0
  if(coop_file_exists(fr_file))then
     nlines = coop_file_numlines(fr_file) 
     if(nlines .gt. 0)then
        call fpsim%open(fr_file, "r")
        do il=1, nlines
           read(fpsim%unit, trim(fmt)) diff
           frmean = frmean+diff
           do j = 0, n
              do i=0, j
                 cov(i, j) = cov(i, j) + diff(i)*diff(j)
              enddo
           enddo
           weight = weight + 1
           if(weight .ge. n_sim) exit
        enddo
        call fpsim%close()
        write(*,*) "Loaded "//trim(coop_num2str(weight))//" lines from checkpoint"
     endif
  endif

  call patch_n%init(stack_type, n, dr, mmax = mmax)
  patch_s = patch_n

  call fpsim%open(fr_file, "a")
  do while(weight .lt. n_sim)
     if(mod(weight, 20) .eq. 0)then
        if(nmaps_wanted .eq. 3)then
           call map%get_fullCls(imask, polmask)
        else
           call map%get_fullCls(imask)
        endif
        sim%cl = map%cl
     endif
     call sim%simulate()
     tmp%map = sim%map
     tmp%ordering = sim%ordering
     weight = weight + 1
     call tmp%get_listpix(listpix, listangle, spot_type, threshold, spots_mask)
     call sim%stack_north_south(patch_n, patch_s, listpix, listangle, hdir, stack_mask)
     write(*,"(A)") "sim #"//trim(coop_num2str(weight))//", num_N = "//trim(coop_num2str(patch_n%nstack_raw))//", num_S = "//trim(coop_num2str(patch_s%nstack_raw))
     call patch_n%get_radial_profile(imap = imap, m = 0)
     call patch_s%get_radial_profile(imap = imap, m = 0)
     diff = patch_n%fr(:, 0, imap) - patch_s%fr(:, 0, imap)
     frmean = frmean + diff
     do j = 0, n
        do i=0, j
           cov(i, j) = cov(i, j) + diff(i) * diff(j)
        enddo
     enddo
     write(fpsim%unit, trim(fmt)) diff
     flush(fpsim%unit)
     write(*, trim(fmtscreen)) diff
  enddo
  call fpsim%close()
  frmean = frmean/weight
  cov = cov/weight
  do j = 0, n
     do i=0, j
        cov(i, j) = cov(i, j) - frmean(i)*frmean(j)
     enddo
  enddo

  do j=0, n
     do i= j+1, n
        cov(i, j) = cov(j, i)
     enddo
     err(j) = sqrt(cov(j,j))
  enddo

  psi = cov
  call coop_matsym_diagonalize(psi, eig, sort = .true.)
  do i=0,n
     if(psi(0,i).lt.0.d0)then
        psi(:,i) = -psi(:,i)
     endif
  enddo


  call fp%open(prefix//"_pca.dat", "w")
  do i=0, n
     write(fp%unit, "("//trim(coop_num2str(ncols_used))//"G16.7)") psi(i, n:n-ncols_used+1:-1)
  enddo
  call fp%close()
  
  do i=0, n
     cov(i,i) =cov(i,i)+epsilon
  enddo

  call coop_matsym_inverse(cov)

  sim = map
  call sim%get_listpix(listpix, listangle, spot_type, threshold, spots_mask)
  call map%stack_north_south(patch_n, patch_s, listpix, listangle, hdir, stack_mask )
  call patch_s%get_radial_profile(imap = imap, m = 0)
  call patch_n%get_radial_profile(imap = imap, m = 0)
  diff = patch_n%fr(:, 0, imap) - patch_s%fr(:, 0, imap)
  chisq =  dot_product(diff, matmul(cov, diff))

  nlines = coop_file_numlines(fr_file) 
  ismall = 0
  ilarge = 0
  call fpsim%open(fr_file, "r")
  do il=1, nlines
     read(fpsim%unit, trim(fmt)) diff
     if (dot_product(diff, matmul(cov, diff)) .ge. chisq)then
        ilarge = ilarge + 1
     else
        ismall = ismall + 1
     endif
  end do
  write(*,*) "probability of seeing chi^2 >= chi^2(data) is ", real(ilarge)/real(nlines)

!!$  do i=n, 1,-1
!!$     print*, sqrt(eig(i)), dot_product(diff, psi(:, i))/sqrt(eig(i))
!!$  enddo
!!$

  call fig%open(prefix//"_radial_profile.txt")
  call fig%init(xlabel = "$r$", ylabel = "$T(\mu K)$")
  call coop_asy_curve(fig, patch_n%r, diff, legend = "NS diff: "//prefix, color="red", linetype = "dashed")
  call coop_asy_curve(fig, patch_n%r, frmean, legend = "NS diff: simu mean", color = "black", linetype = "solid")
  do i = 0, n
     call coop_asy_error_bar(fig, patch_n%r(i), frmean(i), dy_plus = err(i), dy_minus=err(i))
  enddo
  call coop_asy_legend(fig)
  call fig%close()


end program test
