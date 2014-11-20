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
  integer i, n, nmaps, nsims, m_want, map_want, isim, j, start1, end1, start2, end2, n1, n2
  COOP_REAL::dr, pvalue
  type(coop_healpix_maps)::map
  COOP_STRING::prefix
  COOP_REAL, dimension(:,:,:,:),allocatable::f
  COOP_REAL, dimension(:),allocatable::r, rsq, mean, std, chisq
  COOP_REAL, dimension(:,:),allocatable::cov
  COOP_INT::nproj
  COOP_REAL,dimension(:,:),allocatable::vec
  COOP_SHORT_STRING::label
  type(coop_asy)::fig
  
  if(iargc() .lt. 7)then
     write(*,*) "./POSTSST prefix start1 end1 start2 end2 nproj map_want [Label]"
     stop
  endif
  prefix = coop_InputArgs(1)
  start1 = coop_str2int(coop_InputArgs(2))
  end1 = coop_str2int(coop_InputArgs(3))
  start2 = coop_str2int(coop_InputArgs(4))  
  end2 = coop_str2int(coop_InputArgs(5))
  nproj = coop_str2int(coop_InputArgs(6))
  map_want = coop_str2int(coop_InputArgs(7))
  label = trim(coop_InputArgs(8))
  if(trim(label).eq."") label = "P"
  nsims = max(end1, end2)
  n1 = end1- start1+1
  n2 = end2 - start2 + 1
  call fp%open(trim(prefix)//"_info.txt", "r")
  read(fp%unit,*) n, nmaps, dr
  call fp%close()
  nproj = min(n+1, nproj)
  
  dr = dr*coop_SI_arcmin
  allocate(f(0:n, 0:mmax/2, nmaps, 0:nsims), r(0:n), rsq(0:n), mean(nproj), chisq(0:nsims), vec(nproj, 0:nsims), cov(nproj, nproj), std(nproj))
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
  enddo
  f(0:2, 4/2, :, :) = 0.d0  !!remove the m=4 mode generated by pixel effect 
  call fp%close()
  f = f*1.e6

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
        call fr2vec(f(0:n, m_want/2, map_want, i), vec(:, i))
     enddo
     
     do i=1, nproj
        mean(i) = sum(vec(i, start1:end1))/n1
     enddo
     
     do i=1, nproj
        do j=1, i
           cov(i, j) = sum((vec(i, start1:end1)-mean(i))*(vec(j, start1:end1) - mean(j)))/n1
           cov(j, i) = cov(i, j)
        enddo
        cov(i, i) = cov(i, i)*1.0001+1.d-8
        std(i) = sqrt(cov(i,i))
     enddo
     call coop_matsym_inverse(cov)
     do i=0, nsims
        chisq(i) = dot_product(vec(:,i)-mean,  matmul(cov, vec(:, i)-mean))
     enddo
     write(*,"(A, I5, A, F12.4)") "joint; m=", m_want, ", p-value = ", count(chisq(start2:end2).gt.chisq(0))/dble(n2), " chi^2=", chisq(0)
     if(nproj .eq. n + 1)then !!make the scattering plot
        call fig%open(trim(prefix)//"_fig"//trim(label)//COOP_STR_OF(m_want)//".txt")
        call fig%init(xlabel = "$\omega$", ylabel  = "$"//trim(label)//"_"//trim(COOP_STR_OF(m_want))//" (\mu K)$")
        call fig%band(r, mean-std*2, mean+std*2, colorfill = trim(coop_asy_gray_color(0.65)), linecolor = "invisible")
        call fig%band(r, mean-std, mean+std, colorfill = trim(coop_asy_gray_color(0.4)), linecolor = "invisible")        
        call fig%curve(r, f(0:n,m_want/2, map_want, 0), color = "red", linetype = "solid", linewidth = 1.5, legend = "Planck")
        call fig%curve(r, mean, color = "blue", linetype = "dotted", linewidth = 1.5, legend = "FFP8 mean")
        call fig%curve(r, mean+std, color = trim(coop_asy_gray_color(0.4)), linewidth = 0.25, legend = "FFP8 1-$\sigma$")
        call fig%curve(r, mean+std, color = trim(coop_asy_gray_color(0.65)), linewidth = 0.25, legend = "FFP8 2-$\sigma$")        
        call fig%legend(0.7, 0.92)
        call fig%close()
     endif
  end do


contains


  subroutine fr2vec(d, v)
    COOP_REAL d(0:n)
    COOP_REAL v(nproj)
    if(nproj .ge. n+1)then
       v(1:n+1) = d(0:n)
    else
       call coop_chebfit(n+1, rsq, d, nproj, 0.d0, rsq(n), v)
    endif
  end subroutine fr2vec

  
end program test