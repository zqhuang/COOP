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
  COOP_REAL, dimension(:,:,:),allocatable::fr
  COOP_REAL, dimension(:,:),allocatable::f, cov
  COOP_REAL, dimension(:),allocatable::r, mean, std

  type(coop_asy)::fig
  
  if(iargc() .lt. 4)then
     write(*,*) "./POSTHA prefix nsims map_want m_want"
     stop
  endif
  prefix = coop_InputArgs(1)
  nsims = coop_str2int(coop_InputArgs(2))
  map_want = coop_str2int(coop_InputArgs(3))
  m_want =  coop_str2int(coop_InputArgs(4))
  

  call fp%open(trim(prefix)//"_info.txt", "r")
  read(fp%unit,*) n, nmaps, dr
  call fp%close()
  
  dr = dr*coop_SI_arcmin
  allocate(fr(0:n, 0:mmax/2, nmaps), r(0:n), mean(0:n), std(0:n), f(0:n, 0:nsims))
  do i=0,n
     r(i) = (dr*i)
  enddo
  
  call fp%open(trim(prefix)//trim(coop_num2str(i))//".dat", "ru")
  do isim = 0, nsims
     read(fp%unit) j
     if(j.ne.isim)then
        write(*,*) isim, j
        stop "data file broken"
     endif
     read(fp%unit) fr
     f(0:n, isim) = fr(0:n, m_want/2, map_want)
  enddo
  call fp%close()
  do i = 0, n
     mean(i) = sum(f(i, 1:nsims))/nsims
     std(i) = sqrt(sum((f(i, 1:nsims) - mean(i))**2)/nsims)
  enddo
  
  call fig%open(trim(prefix)//"_fig.txt")
  call fig%init(xlabel = "$r$", ylabel  = "radial profile")
  call fig%band(r, mean-std*2, mean+std*2, colorfill = "lightgray", linecolor = "invisible")
  call fig%band(r, mean-std, mean+std, colorfill = "gray", linecolor = "invisible")
  call fig%curve(r, f(:,0), color = "red", linetype = "solid", linewidth = 1.5, legend = "data")
  call fig%curve(r, mean, color = "blue", linetype = "dotted", linewidth = 1.5, legend = "sim. mean")
  
  call fig%legend(0.1, 0.95, 1)
  call fig%close()
  
end program test
