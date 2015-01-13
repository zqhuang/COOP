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
  type(coop_file)::fp
  COOP_INT,parameter::mmax = 4
  COOP_REAL,parameter::normK = 1.d6
  integer i, n, nmaps, nsims, m_want, map_want, isim, j, start1, end1, start2, end2, n1, n2
  COOP_REAL::dr, pvalue
  type(coop_healpix_maps)::map
  COOP_STRING::prefix
  COOP_REAL, dimension(:,:,:,:),allocatable::f
  COOP_REAL, dimension(:,:,:), allocatable:: fsim_base, fdata_base
  COOP_REAL, dimension(:),allocatable::r, rsq, mean, chisq, tmpf
  COOP_REAL, dimension(:,:),allocatable::cov, bounds
  COOP_INT::nproj
  COOP_REAL::chisq0
  COOP_REAL,dimension(:,:),allocatable::vec
  COOP_SHORT_STRING::label
  COOP_STRING::bounds_file
  COOP_STRING::data_mean_file = "planck14/data_mean_file.dat"
  COOP_STRING::sim_mean_file = "ffp8/sim_mean_file.dat"  
  type(coop_asy)::fig
  if(iargc() .lt. 9)then
     write(*,*) "./POSTSST data_mean_file, sim_mean_file, prefix start1 end1 start2 end2 nproj map_want [Label] [Bounds_File]"
     stop
  endif
  data_mean_file = trim(coop_InputArgs(1))
  sim_mean_file = trim(coop_InputArgs(2))  
  prefix = coop_InputArgs(3)
  start1 = coop_str2int(coop_InputArgs(4))
  end1 = coop_str2int(coop_InputArgs(5))
  start2 = coop_str2int(coop_InputArgs(6))  
  end2 = coop_str2int(coop_InputArgs(7))
  nproj = coop_str2int(coop_InputArgs(8))
  map_want = coop_str2int(coop_InputArgs(9))
  label = trim(coop_InputArgs(10))
  bounds_file = trim(coop_InputArgs(11))
  
  if(trim(label).eq."") label = "P"

  nsims = max(end1, end2)
  n1 = end1- start1+1
  n2 = end2 - start2 + 1
  
  call fp%open(trim(prefix)//"_info.txt", "r")
  read(fp%unit,*) n, nmaps, dr
  call fp%close()

  
  
  nproj = min(n+1, nproj)
  
  dr = dr*coop_SI_arcmin
  allocate(f(0:n, 0:mmax/2, nmaps, 0:nsims), fsim_base(0:n, 0:mmax/2, nmaps),  fdata_base(0:n, 0:mmax/2, nmaps), r(0:n), rsq(0:n), mean(nproj), chisq(0:nsims), vec(nproj, 0:nsims), cov(nproj, nproj), bounds(-2:2, nproj), tmpf(n) )
  call fp%open(sim_mean_file, "r")
  read(fp%unit, *) fsim_base
  call fp%close()

  call fp%open(data_mean_file, "r")
  read(fp%unit, *) fdata_base
  call fp%close()
  
  do i=0,n
     r(i) = (dr*i)
  enddo
  rsq = r**2
  call fp%open(trim(prefix)//".dat", "ru")
  do isim = 0, nsims
     read(fp%unit) j
     if(j.ne.isim)then
        write(*,*) isim, j
        stop "data file broken"
     endif
     read(fp%unit) f(0:n, 0:mmax/2, 1:nmaps, isim)
     if(isim.eq.0)then
        f(0:n, 0:mmax/2, 1:nmaps, isim) = f(0:n, 0:mmax/2, 1:nmaps, isim)*normK - fsim_base        
     else
        f(0:n, 0:mmax/2, 1:nmaps, isim) = f(0:n, 0:mmax/2, 1:nmaps, isim)*normK - fdata_base
     endif
  enddo
  call fp%close()
  call fp%close()
  

  !!collapse to complex radial modes
  if(map_want .gt. nmaps)then
     do isim = 0, nsims
        do m_want = 0, mmax, 2
           do i = 0, n
              f(i, m_want/2, 1, isim) = sum(f(i, m_want/2, 1:nmaps, isim))/nmaps
           enddo
        enddo
     enddo
     map_want = 1     
  endif
  do m_want = 0, mmax, 2
     do i = 0, nsims
        call fr2vec(f(0:n, m_want/2, map_want, i), vec(:, i), m_want)
     enddo
     do i=1, nproj
        mean(i) = sum(vec(i, start1:end1))/n1
     enddo
     
     if(trim(bounds_file).ne."")then
        call fp%open(trim(bounds_file)//COOP_STR_OF(m_want)//".txt", "r")
        read(fp%unit, *) bounds
        call fp%close()
     else
        do i=1, nproj
           call coop_get_bounds( vec(i, start1:end1), (/ 0.023d0, 0.1585d0, 0.5d0, 0.8415d0, 0.977d0 /), bounds(-2:2, i) )
        enddo
     endif
     
     do i=1, nproj
        do j=1, i
           cov(i, j) = sum((vec(i, start1:end1)-mean(i))*(vec(j, start1:end1) - mean(j)))/n1
           cov(j, i) = cov(i, j)
        enddo
        cov(i, i) = cov(i, i)
     enddo
     call coop_matsym_inverse(cov)
     if(nproj .le. 5)then
        write(*,*) "mean:"
        write(*,"(5G14.5)") mean
        write(*,*) "cov:"        
        do i=1, nproj
           write(*,"(5G14.5)") cov(i, :)
        enddo
     endif
     do i=0, nsims
        chisq(i) = dot_product(vec(:,i)-mean,  matmul(cov, vec(:, i)-mean))
     enddo
     write(*,"(A, I5, A, F12.4, A, F12.4)") "m=", m_want, ", p-value = ", count(chisq(start2:end2).gt.chisq(0))/dble(n2), ", chi^2=", chisq(0)
     chisq0 = dot_product(vec(:,0),  matmul(cov, vec(:, 0)))
     write(*,"(A, I5, A, F12.4, A, F12.4)") "    Noise-free case: p-value = ", count(chisq(start2:end2).gt.chisq0)/dble(n2), ", chi^2=", chisq0
     
     if(nproj .eq. n + 1)then !!make the scattering plot
        call fig%open(trim(prefix)//"_fig"//trim(label)//COOP_STR_OF(m_want)//".txt")
        call fig%init(xlabel = "$\varpi$", ylabel  = "$\delta "//trim(label)//"_"//trim(COOP_STR_OF(m_want))//" (\mu K)$")
        call fig%band(r, bounds(-2,:)-bounds(0,:), bounds(2,:)-bounds(0,:), colorfill = trim(coop_asy_gray_color(0.65)), linecolor = "invisible")
        call fig%band(r, bounds(-1,:)-bounds(0,:), bounds(1,:)-bounds(0,:), colorfill = trim(coop_asy_gray_color(0.42)), linecolor = "invisible")        
        do i=1, 20
           j = coop_random_index(nsims)
           if(i.eq.1)then
              call fig%curve(r, vec(:, j) - bounds(0, :), color = "HEX:1010EF", linetype = "dotted", linewidth = 1., legend = "FFP8 sims.")
           else
              call fig%curve(r, vec(:, j) - bounds(0, :), color = "HEX:1010EF", linetype = "dotted", linewidth = 1.)
           endif
        enddo
        call fig%curve(r, vec(:, 0) - bounds(0, :), color = "HEX:21C0FC", linetype = "solid", linewidth = 2., legend = "Planck")        
        call fig%add_legend(legend = "1-$\sigma$ bound", color = trim(coop_asy_gray_color(0.42)))
        call fig%add_legend(legend = "2-$\sigma$ bound", color = trim(coop_asy_gray_color(0.65)))        
        call fig%expand(0., 0., 0.01, 0.25)
        call fig%legend(0.05, 0.95, 2)
        call fig%close()

        call fig%open(trim(prefix)//"_full"//trim(label)//COOP_STR_OF(m_want)//".txt")
        call fig%init(xlabel = "$\varpi$", ylabel  = "$ "//trim(label)//"_"//trim(COOP_STR_OF(m_want))//" (\mu K)$")
        call fig%band(r, bounds(-2,:), bounds(2,:), colorfill = trim(coop_asy_gray_color(0.65)), linecolor = "invisible")
        call fig%band(r, bounds(-1,:), bounds(1,:), colorfill = trim(coop_asy_gray_color(0.42)), linecolor = "invisible")        
        call fig%curve(r, vec(:, 0), color = "HEX:21C0FC", linetype = "solid", linewidth = 1.8, legend = "Planck")
        call fig%curve(r, bounds(0,:), color = "red", linetype = "dashed", linewidth = 1.8, legend = "mean")                        
        call fig%add_legend(legend = "1-$\sigma$ bound", color = trim(coop_asy_gray_color(0.42)))
        call fig%add_legend(legend = "2-$\sigma$ bound", color = trim(coop_asy_gray_color(0.65)))        
        call fig%expand(0., 0., 0.01, 0.25)
        call fig%legend(0.05, 0.95, 2)
        call fig%close()
        
     endif
     if(trim(bounds_file).eq."")then
        call fp%open(trim(prefix)//"_median"//trim(label)//COOP_STR_OF(m_want)//".txt", "w")
        write(fp%unit, "("//COOP_STR_OF(size(bounds))//"E16.7)") bounds
        call fp%close()
     endif
  end do


contains

  subroutine fr2vec(d, v, m) !! this is doing simple binning
    COOP_REAL d(0:n)
    COOP_REAL v(nproj), w, sumw
    COOP_INT m, i, n_per_bin, j
    n_per_bin = (n+1)/nproj
    v = 0.d0
    do i = 1, nproj
       sumw = 0.d0
       do j = n_per_bin*(i-1), n_per_bin*i-1
          w = max(0.25d0, dble(j))
          v(i) = v(i) + w*d(j)
          sumw = sumw + w
       enddo
       v(i) = v(i)/sumw
    enddo
  end subroutine fr2vec


  subroutine fr2vec2(d, v, m)  !!this is doing chebyshev fitting
    COOP_REAL d(0:n)
    COOP_REAL v(nproj)
    COOP_INT m
    if(nproj .ge. n+1)then
       v(1:n+1) = d(0:n)
    else
       if(m.eq.0)then
          call coop_chebfit(n+1, rsq, d, nproj, 0.d0, rsq(n), v)
          return
       endif
       tmpf = d(1:n)/r(1:n)**m
       call coop_chebfit(n, rsq(1:n), tmpf, nproj, 0.d0, rsq(n), v)       
    endif
  end subroutine fr2vec2

  
end program test
