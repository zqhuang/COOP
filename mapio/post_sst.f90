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
  integer i, n, nmaps, nsims, m_want, map_want, isim, j
  COOP_REAL::dr
  type(coop_healpix_maps)::map
  COOP_STRING::prefix
  COOP_REAL, dimension(:,:,:,:),allocatable::f
  COOP_REAL, dimension(:),allocatable::r, mean, std, chisq
  COOP_REAL, dimension(:,:),allocatable::cov

  type(coop_asy)::fig
  
  if(iargc() .lt. 2)then
     write(*,*) "./POSTSST prefix nsims label_map1, label_map2 ..."
     stop
  endif
  prefix = coop_InputArgs(1)
  nsims = coop_str2int(coop_InputArgs(2))
  

  call fp%open(trim(prefix)//"_info.txt", "r")
  read(fp%unit,*) n, nmaps, dr
  call fp%close()
  
  if(iargc() .lt. 2 + nmaps)then
     write(*,*) "./POSTSST prefix nsims label_map1, label_map2 ..."
     stop
  endif
  
  dr = dr*coop_SI_arcmin
  allocate(f(0:n, 0:mmax/2, nmaps, 0:nsims), r(0:n), mean(0:n), std(0:n), cov(0:n, 0:n), chisq(0:nsims))
  do i=0,n
     r(i) = (dr*i)
  enddo
  
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
  do map_want = 1, nmaps
     do m_want = 0, mmax, 2
        do i = 0, n
           mean(i) = sum(f(i, m_want/2, map_want, 1:nsims))/nsims
           std(i) = sqrt(sum((f(i, m_want/2, map_want, 1:nsims) - mean(i))**2)/nsims)
        enddo
        do i=0, n
           do j = 0, i
              cov(i, j) = sum((f(i, m_want, map_want, 1:nsims)-mean(i))*(f(j, m_want, map_want, 1:nsims)-mean(j)))/nsims
              cov(j, i) = cov(i, j)
           enddo
        enddo
        call coop_matsym_inverse(cov)
        do i=0, nsims
           chisq(i) = dot_product(f(:, m_want, map_want, i)-mean,  matmul(cov, f(:, m_want, map_want, i)-mean))
        enddo
        write(*,"(A,I5, A, I5, A, F12.4)") "map#",map_want, "; m=", m_want, ", rareness = ", count(chisq(1:nsims).gt.chisq(0))/nsims
        
        call fig%open(trim(prefix)//"_fig"//COOP_STR_OF(map_want)//COOP_STR_OF(m_want)//".txt")
        call fig%init(xlabel = "$\omega$", ylabel  = "$"//trim(coop_InputArgs(2+map_want))//" (\mu K)$")
        call fig%band(r, mean-std*2, mean+std*2, colorfill = trim(coop_asy_gray_color(0.65)), linecolor = "invisible")
        call fig%band(r, mean-std, mean+std, colorfill = trim(coop_asy_gray_color(0.4)), linecolor = "invisible")        
        call fig%curve(r, f(0:n,m_want/2, map_want, 0), color = "red", linetype = "solid", linewidth = 1.5, legend = "Planck")
        call fig%curve(r, mean, color = "blue", linetype = "dotted", linewidth = 1.5, legend = "FFP8 mean")
        call fig%curve(r, mean+std, color = trim(coop_asy_gray_color(0.4)), linewidth = 0.25, legend = "FFP8 1-$\sigma$")
        call fig%curve(r, mean+std, color = trim(coop_asy_gray_color(0.65)), linewidth = 0.25, legend = "FFP8 2-$\sigma$")        
        if(m_want.le.2)then
           call fig%label("$m = "//COOP_STR_OF(m_want)//"$ mode", 0.3, 0.6)
        else
           call fig%label("$m = "//COOP_STR_OF(m_want)//"$ mode", 0.3, 0.4)           
        endif
        call fig%legend(0.7, 0.92)
        call fig%close()
     end do
  enddo

  !!collapse to complex radial modes
  if(nmaps .eq. 2)then
     do isim = 0, nsims
        do m_want = 0, mmax, 2
           f(:, m_want/2, 1, isim) = (f(:, m_want/2, 1, isim)+f(:, m_want/2, 2, isim))/2.d0
        enddo
     enddo
     map_want = 1
     do m_want = 0, mmax, 2
        do i = 0, n
           mean(i) = sum(f(i, m_want/2, map_want, 1:nsims))/nsims
           std(i) = sqrt(sum((f(i, m_want/2, map_want, 1:nsims) - mean(i))**2)/nsims)
        enddo
        do i=0, n
           do j = 0, i
              cov(i, j) = sum((f(i, m_want, map_want, 1:nsims)-mean(i))*(f(j, m_want, map_want, 1:nsims)-mean(j)))/nsims
              cov(j, i) = cov(i, j)
           enddo
        enddo
        call coop_matsym_inverse(cov)
        do i=0, nsims
           chisq(i) = dot_product(f(:, m_want, map_want, i)-mean,  matmul(cov, f(:, m_want, map_want, i)-mean))
        enddo
        write(*,"(A, I5, A, F12.4)") "joint; m=", m_want, ", rareness = ", count(chisq(1:nsims).gt.chisq(0))/nsims
        call fig%open(trim(prefix)//"_figp"//COOP_STR_OF(m_want)//".txt")
        call fig%init(xlabel = "$\omega$", ylabel  = "$ P_"//trim(COOP_STR_OF(m_want))//" (\mu K)$")
        call fig%band(r, mean-std*2, mean+std*2, colorfill = trim(coop_asy_gray_color(0.65)), linecolor = "invisible")
        call fig%band(r, mean-std, mean+std, colorfill = trim(coop_asy_gray_color(0.4)), linecolor = "invisible")        
        call fig%curve(r, f(0:n,m_want/2, map_want, 0), color = "red", linetype = "solid", linewidth = 1.5, legend = "Planck")
        call fig%curve(r, mean, color = "blue", linetype = "dotted", linewidth = 1.5, legend = "FFP8 mean")
        call fig%curve(r, mean+std, color = trim(coop_asy_gray_color(0.4)), linewidth = 0.25, legend = "FFP8 1-$\sigma$")
        call fig%curve(r, mean+std, color = trim(coop_asy_gray_color(0.65)), linewidth = 0.25, legend = "FFP8 2-$\sigma$")        
        if(m_want.le.2)then
           call fig%label("$m = "//COOP_STR_OF(m_want)//"$ mode", 0.3, 0.6)
        else
           call fig%label("$m = "//COOP_STR_OF(m_want)//"$ mode", 0.3, 0.4)           
        endif
        call fig%legend(0.7, 0.92)
        call fig%close()
     end do
  endif
end program test
