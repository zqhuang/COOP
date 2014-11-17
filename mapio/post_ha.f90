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
  COOP_INT::m_want = 0
  integer i, nside, id, pix, iminprob, n, nmaps, nsims, isim, j, ncut, i1, i2, cgt, pix_want
  type(coop_healpix_maps)::map
  COOP_STRING::prefix
  COOP_REAL theta, phi, chisq(0:95), prob(0:95), minprob, l, b, dr
  COOP_REAL,dimension(:,:),allocatable::cov, vec, diff, diffmin, vecmin
  COOP_REAL, dimension(:,:,:),allocatable::frn, frs
  COOP_REAL, dimension(:),allocatable:: rsq, r, mean, fitdiff
  COOP_REAL xlab, ylab
  type(coop_asy)::fig
  logical::single_pix = .false.
  
  if(iargc().lt.4)then
     write(*,*) "./POSTHA prefix nsims ncut m_want"
     write(*,*) " or "     
     write(*,*) "./POSTHA prefix 0 pix_want ncut m_want"
     stop
  endif
  prefix = coop_InputArgs(1)  
  nsims = coop_str2int(coop_InputArgs(2))
  if(nsims .eq. 0)then
     single_pix = .true.
     pix_want = coop_str2int(coop_InputArgs(3))
     ncut = coop_str2int(coop_InputArgs(4))
     m_want = coop_str2int(coop_InputArgs(5))     
  else
     single_pix = .false.
     ncut = coop_str2int(coop_InputArgs(3))
     m_want = coop_str2int(coop_InputArgs(4))
  endif
  call fp%open(trim(prefix)//"info.txt", "r")
  read(fp%unit,*) n, nmaps, dr
  call fp%close()
  dr = dr*coop_SI_arcmin
  allocate(frn(0:n, 0:mmax/2, nmaps), frs(0:n, 0:mmax/2, nmaps),diff(0:n, 0:nsims),diffmin(0:n, 0:nsims), r(0:n), rsq(0:n), vec(ncut, 0:nsims), vecmin(ncut, 0:nsims), cov(ncut, ncut), mean(ncut), fitdiff(0:n))
  do i=0,n
     r(i) = (dr*i)
  enddo
  rsq = r**2
  
  call map%init(nside = 4, nmaps=1, spin = (/ 0 /))
  map%map = 0.
  minprob = 10.
  iminprob = 0
  if(single_pix)then
     i = pix_want
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
        if(nmaps .eq. 1)then
           diff(0:n, isim) = (frn(0:n, m_want/2, 1)- frs(0:n, m_want/2, 1))*1.d6
        elseif(nmaps.eq.2)then
           diff(0:n, isim) = (frn(0:n, m_want/2, 1) + frn(0:n, m_want/2, 2) - frs(0:n, m_want/2, 1) - frs(0:n, m_want/2, 2) )*0.5d6
        else
           stop "unknown nmaps"
        endif
        call fr2vec(diff(:, isim),vec(:, isim))
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
     minprob = prob(i)
     iminprob = i
     diffmin = diff
     vecmin = vec
     call pix2ang_ring(map%nside, i, theta, phi)
     theta = coop_pi - theta
     phi = coop_pi + phi
     call ang2pix_ring(map%nside, theta, phi, pix)
     map%map(pix, 1) = map%map(i, 1)
     print*, i, prob(i), chisq(i)
  else
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
           if(nmaps .eq. 1)then
              diff(0:n, isim) = (frn(0:n, m_want/2, 1)- frs(0:n, m_want/2, 1))*1.d6
           elseif(nmaps.eq.2)then
              diff(0:n, isim) = (frn(0:n, m_want/2, 1) + frn(0:n, m_want/2, 2) - frs(0:n, m_want/2, 1) - frs(0:n, m_want/2, 2) )*0.5d6
           else
              stop "unknown nmaps"
           endif
           call fr2vec(diff(:, isim),vec(:, isim))
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
           diffmin = diff
           vecmin = vec
        endif
        call pix2ang_ring(map%nside, i, theta, phi)
        theta = coop_pi - theta
        phi = coop_pi + phi
        call ang2pix_ring(map%nside, theta, phi, pix)
        map%map(pix, 1) = map%map(i, 1)
        print*, i, prob(i), chisq(i)
     enddo
  endif
  call map%write(trim(prefix)//"powercut"//trim(coop_num2str(ncut))//"_probs_m"//COOP_STR_OF(m_want)//".fits")
  call pix2ang_ring(map%nside, iminprob, theta, phi)
  call coop_healpix_ang2lb(theta, phi, l, b)
  write(*,*) "min prob = ", minprob
  if(b.gt.0.d0)then
     l = l + 180.d0
     b =  - b
  endif
  write(*,*) "direction l = ", nint(l), " b = ", nint(b)
  call fig%open(trim(prefix)//"powercut"//trim(coop_num2str(ncut))//"_fr_m"//COOP_STR_OF(m_want)//".txt")
  call fig%init(xlabel = "$\omega$", ylabel  = "$\delta f (\mu K)$")
  call fig%curve(r, diffmin(:, 0), color = "red", linetype = "solid", linewidth = 2., legend = "Planck")
  i = 1
  call fig%curve(r, diffmin(:, i), color = "gray", linetype = "dashed", linewidth = 0.5, legend = "FFP8")
  do i = 2, nsims, nsims/25
     call fig%curve(r, diffmin(:, i), color = "gray", linetype = "dashed", linewidth = 0.5)
  enddo
  do i = 0, n
     call coop_chebeval(ncut, 0.d0, rsq(n), vecmin(:, 0), rsq(i), fitdiff(i))
  enddo
  call fig%curve(r, fitdiff, color = "blue", linetype = "dotted", linewidth = 1.5, legend = "fit")
  call fig%legend(0.1, 0.95, 1)
  call fig%close()
  call fig%open(trim(prefix)//"powercut"//COOP_STR_OF(ncut)//"_distr_m"//COOP_STR_OF(m_want)//".txt")
  call fig%init(xlabel = "$c_0$", ylabel = "$c_1$")
  call fig%dots(vecmin(1, 1:nsims), vecmin(2, 1:nsims), "black")
  call fig%dot(vecmin(1, 0), vecmin(2, 0), "red", "x")
  ylab = maxval(vecmin(2, 0:nsims))
  xlab = maxval(vecmin(1, 0:nsims))
  call fig%dot(xlab, ylab*1.1, "invisible", "x")
  call coop_asy_label(fig, "x   Planck", xlab*0.2, ylab*1.08, "red")
  call coop_asy_label(fig, "$\bullet$   FFP8", xlab*0.2, ylab*1.02, "black")
  call fig%close()
contains

  subroutine fr2vec(d, v)
    COOP_REAL d(0:n)
    COOP_REAL v(ncut)
    call coop_chebfit(n+1, rsq, d, ncut, 0.d0, rsq(n), v)
  end subroutine fr2vec
  
end program test
