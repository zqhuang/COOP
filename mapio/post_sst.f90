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
  COOP_REAL, dimension(:),allocatable::r, mean, std

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
  allocate(f(0:n, 0:mmax/2, nmaps, 0:nsims), r(0:n), mean(0:n), std(0:n))
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
  call fp%close()
  f = f*1.e6
  do map_want = 1, nmaps
     do m_want = 0, mmax, 2
        do i = 0, n
           mean(i) = sum(f(i, m_want/2, map_want, 1:nsims))/nsims
           std(i) = sqrt(sum((f(i, m_want/2, map_want, 1:nsims) - mean(i))**2)/nsims)
        enddo
        
        call fig%open(trim(prefix)//"_fig"//COOP_STR_OF(map_want)//COOP_STR_OF(m_want)//".txt")
        call fig%init(xlabel = "$r$", ylabel  = "$"//trim(coop_InputArgs(2+map_want))//" (\mu K)$")
        call fig%band(r, mean-std*2, mean+std*2, colorfill = trim(coop_asy_gray_color(0.75)), linecolor = "invisible")
        call fig%curve(r, mean+std, color = trim(coop_asy_gray_color(0.75)), linewidth = 0.5, legend = "FFP8 2-$\sigma$")
        call fig%band(r, mean-std, mean+std, colorfill = "gray", linecolor = "invisible")
        call fig%curve(r, mean+std, color = "gray", linewidth = 0.5, legend = "FFP8 1-$\sigma$")
        call fig%curve(r, f(:,m_want, map_want, 0), color = "red", linetype = "solid", linewidth = 1.5, legend = "Data")
        call fig%curve(r, mean, color = "blue", linetype = "dotted", linewidth = 1.5, legend = "FFP8 mean")
        call fig%label("$m = "//COOP_STR_OF(m_want)//"$ mode", 0.4, 0.6)
        call fig%legend(0.6, 0.92)
        call fig%close()
     end do
  enddo
  
end program test
