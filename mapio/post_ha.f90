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
  COOP_INT,parameter::m_want = 0
  COOP_INT,parameter::map_want = 1
  integer i, nside, id, pix, iminprob, n, nmaps, nsims, isim, j, ncut, i1, i2, cgt
  type(coop_healpix_maps)::map
  COOP_STRING::prefix
  COOP_REAL theta, phi, chisq(0:95), prob(0:95), minprob, l, b, dr
  COOP_REAL,dimension(:,:),allocatable::cov, vec 
  COOP_REAL, dimension(:,:,:),allocatable::frn, frs
  COOP_REAL, dimension(:),allocatable::diff, rsq, r, mean, datadiff

  type(coop_asy)::fig
  
  prefix = coop_InputArgs(1)
  if(trim(prefix).eq."")then
     write(*,*) "./POSTHA prefix nsims ncut"
     stop
  endif
  nsims = coop_str2int(coop_InputArgs(2))
  ncut = coop_str2int(coop_InputArgs(3))

  call fp%open(trim(prefix)//"info.txt", "r")
  read(fp%unit,*) n, nmaps, dr
  call fp%close()
  dr = dr*coop_SI_arcmin
  allocate(frn(0:n, 0:mmax/2, nmaps), frs(0:n, 0:mmax/2, nmaps),diff(0:n), datadiff(0:n), r(0:n), rsq(0:n), vec(ncut, 0:nsims), cov(ncut, ncut), mean(ncut))
  do i=0,n
     r(i) = (dr*i)
  enddo
  rsq = r**2
  
  call map%init(nside = 4, nmaps=1, spin = (/ 0 /))
  map%map = 0.
  minprob = 10.
  iminprob = 0
  do i=0, 95
     write(*,*) "reading "//trim(prefix)//trim(coop_num2str(i))//".dat"
     call fp%open(trim(prefix)//trim(coop_num2str(i))//".dat", "ru")
     do isim = 0, nsims
        read(fp%unit) j
        if(j.ne.isim)then
           write(*,*) isim, j
           stop "data file broken"
        endif
        read(fp%unit) frn
        read(fp%unit) frs
        if(isim.eq.0)   datadiff = frn(0:n, m_want/2, map_want)- frs(0:n, m_want/2, map_want)
        call fr2vec(vec(:, isim))
     enddo
     call fp%close()
     cov = 0.d0
     do i1 = 1, ncut
        mean(i1) = sum(vec(i1, 1:nsims))/nsims
     enddo
     do i1 = 1, ncut     
        do i2 = 1, i1
           cov(i1,i2) = sum((vec(i1, 1:nsims)-mean(i1))*(vec(i2, 1:nsims)-mean(i2)))/nsims
           cov(i2,i1) = cov(i1,i2)
        enddo
     enddo
     call coop_matsym_inverse(cov)     
     chisq(i) = dot_product(vec(:, 0) - mean, matmul(cov, vec(:, 0) - mean))
     cgt = 0
     do isim = 1, nsims
        if(dot_product(vec(:, isim) - mean, matmul(cov, vec(:, isim) - mean)).gt. chisq(i)) cgt = cgt + 1
     enddo
     prob(i) = dble(cgt) / nsims
     map%map(i, 1) = log10(max(prob(i), 1.d-4))
     if(minprob .gt. prob(i))then
        minprob = prob(i)
        iminprob = i
     endif
     call pix2ang_ring(map%nside, i, theta, phi)
     theta = coop_pi - theta
     phi = coop_pi + phi
     call ang2pix_ring(map%nside, theta, phi, pix)
     map%map(pix, 1) = map%map(i, 1)
     print*, i, prob(i), chisq(i)
  enddo
  call map%write(trim(prefix)//"_probs.fits")
  call pix2ang_ring(map%nside, iminprob, theta, phi)
  call coop_healpix_ang2lb(theta, phi, l, b)
  write(*,*) "min prob = ", minprob
  if(b.gt.0.d0)then
     l = l + 180.d0
     b =  - b
  endif
  write(*,*) "direction l = ", nint(l), " b = ", nint(b)
!!$  call coop_asy_histogram(chisq, 10, "chisq_hist.txt")
!!$  call coop_asy_histogram(log(max(prob, 1.d-4)), 10, "logprob_hist.txt")
  call fig%open(trim(prefix)//"_profile_fit.txt")
  call fig%init(xlabel = "$r$", ylabel  = "\delta T(\mu K)")
  call coop_asy_curve(fig, r, datadiff, color = "red", linetype = "solid", linewidth = 1.5)
  do i = 1, ncut
     call coop_chebeval(ncut, 0.d0, rsq(n), vec(:, 0), rsq(i), diff(i))
  enddo
  call fig%curve(r, diff, color = "blue", linetype = "dotted", linewidth = 1.3)
  call fig%legend(0.1, 0.1, 1)
  call fig%close()
contains

  subroutine fr2vec(v)
    COOP_REAL v(ncut)
    diff = frn(0:n, m_want/2, map_want)- frs(0:n, m_want/2, map_want)
    call coop_chebfit(n+1, rsq, diff, ncut, 0.d0, rsq(n), v)
  end subroutine fr2vec
  
end program test
