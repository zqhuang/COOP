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
  integer i, nside, id, pix, iminprob, n, nmaps, nsims, isim, j, ncut, i1, i2, cgt, pix_want, ind_data
  type(coop_healpix_maps)::map
  COOP_STRING::prefix
  COOP_REAL theta, phi, chisq(0:95), prob(0:95), minprob, l, b, dr
  COOP_INT::start1, end1, start2, end2, n1, n2
  COOP_REAL,dimension(:,:),allocatable::cov, vec, diff, diffmin, vecmin
  COOP_REAL, dimension(:,:,:),allocatable::frn, frs
  COOP_REAL, dimension(:),allocatable:: rsq, r, mean, fitdiff
  type(coop_asy)::fig
  logical::single_pix = .false.
  call coop_MPI_init()
  call coop_random_init()
  if(iargc().lt.7)then
     write(*,*) "./POSTHA prefix start end1 start2 end2 ncut m_want"
     write(*,*) " or "     
     write(*,*) "./POSTHA prefix start1, end1, start2, end2 ncut m_want pix_want"
     stop
  endif
  prefix = coop_InputArgs(1)  
  start1 = coop_str2int(coop_InputArgs(2))
  end1 =  coop_str2int(coop_InputArgs(3))
  start2 = coop_str2int(coop_InputArgs(4))
  end2 = coop_str2int(coop_InputArgs(5))
  ncut = coop_str2int(coop_InputArgs(6))
  m_want = coop_str2int(coop_InputArgs(7))
  nsims = max(end1, end2)
  n1 = end1 - start1 + 1 
  n2 = end2 - start2 + 1
  if(iargc() .ge. 8)then
     single_pix = .true.
     pix_want = coop_str2int(coop_InputArgs(8))
  else
     single_pix = .false.
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
        mean(i1) = sum(vec(i1, start1:end1))/n1
     enddo
     do i1 = 1, ncut     
        do i2 = 1, i1
           cov(i1,i2) = sum((vec(i1, start1:end1)-mean(i1))*(vec(i2, start1:end1)-mean(i2)))/n1
           cov(i2,i1) = cov(i1,i2)
        enddo
     enddo
     
     call coop_matsym_inverse(cov)     
     chisq(i) = dot_product(vec(:, 0) - mean, matmul(cov, vec(:, 0) - mean))
     cgt = 0
     do isim = start2, end2
        if(dot_product(vec(:, isim) - mean, matmul(cov, vec(:, isim) - mean)).gt. chisq(i)) cgt = cgt + 1
     enddo
     prob(i) = dble(cgt) / n2
     map%map(i, 1) = log10(max(prob(i), 1.d-4))
     minprob = prob(i)
     iminprob = i
     diffmin = diff
     vecmin = vec
     call pix2ang_ring(map%nside, i, theta, phi)
     call coop_healpix_ang2lb(theta, phi, l, b)
     if(b.gt.0.d0)then
        l = l + 180.d0
        b =  - b
     endif
     theta = coop_pi - theta
     phi = coop_pi + phi
     call ang2pix_ring(map%nside, theta, phi, pix)
     map%map(pix, 1) = map%map(i, 1)
     write(*, "(3I5, 2F10.4)")i, nint(l), nint(b), prob(i), chisq(i)
  else
     ind_data = 0
     do i=0, 95
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
           mean(i1) = sum(vec(i1, start1:end1))/n1
        enddo
        do i1 = 1, ncut     
           do i2 = 1, i1
              cov(i1,i2) = sum((vec(i1, start1:end1)-mean(i1))*(vec(i2, start1:end1)-mean(i2)))/n1
              cov(i2,i1) = cov(i1,i2)
           enddo
        enddo
        call coop_matsym_inverse(cov)     
        chisq(i) = dot_product(vec(:, 0) - mean, matmul(cov, vec(:, 0) - mean))
        cgt = 0
        do isim = start2, end2
           if(dot_product(vec(:, isim) - mean, matmul(cov, vec(:, isim) - mean)).gt. chisq(i)) cgt = cgt + 1
        enddo
        prob(i) = dble(cgt) / n2
        map%map(i, 1) = log10(max(prob(i), 1.d-4))
        if(minprob .gt. prob(i))then
           minprob = prob(i)
           iminprob = i
           diffmin = diff
           vecmin = vec
        endif
        call pix2ang_ring(map%nside, i, theta, phi)
        call coop_healpix_ang2lb(theta, phi, l, b)
        if(b.gt.0.d0)then
           l = l + 180.d0
           b =  - b
        endif        
        theta = coop_pi - theta
        phi = coop_pi + phi
        call ang2pix_ring(map%nside, theta, phi, pix)
        map%map(pix, 1) = map%map(i, 1)
        write(*, "(3I5, 2F10.4)")i, nint(l), nint(b), prob(i), chisq(i)
     enddo
  endif
  call map%write(trim(prefix)//"powercut"//trim(coop_num2str(ncut))//"_"//COOP_STR_OF(nsims)//"sims_probs_m"//COOP_STR_OF(m_want)//".fits")
  call pix2ang_ring(map%nside, iminprob, theta, phi)
  call coop_healpix_ang2lb(theta, phi, l, b)
  if(b.gt.0.d0)then
     l = l + 180.d0
     b =  - b
  endif
  write(*,*) "prob = ", minprob  
  write(*,*) "direction l = ", nint(l), " b = ", nint(b), "ipix = ", iminprob
  call fig%open(trim(prefix)//"powercut"//trim(coop_num2str(ncut))//"_"//COOP_STR_OF(nsims)//"sims_fr_m"//COOP_STR_OF(m_want)//".txt")
  call fig%init(xlabel = "$\varpi$", ylabel  = "$\delta f (\mu K)$")
  call fig%curve(r, diffmin(:, 0), color = "HEX:21C0FC", linetype = "solid", linewidth = 2., legend = "Planck")
  do i = 0, n
     call coop_chebeval(ncut, 0.d0, rsq(n), vecmin(:, 0), rsq(i), fitdiff(i))
  enddo
  if(m_want.gt.0) fitdiff = fitdiff*r**m_want
  call fig%curve(r, fitdiff, color = "blue", linetype = "dotted", linewidth = 1.5, legend = "Planck fit")
  i = 1
  call fig%curve(r, diffmin(:, i), color =trim(coop_asy_gray_color(0.2)), linetype = "dashed", linewidth = 0.5, legend = "FFP8")
  do j = 0, n
     call coop_chebeval(ncut, 0.d0, rsq(n), vecmin(:, i), rsq(j), fitdiff(j))
  enddo
  if(m_want.gt.0) fitdiff = fitdiff*r**m_want  
  call fig%curve(r, fitdiff, color = trim(coop_asy_gray_color(0.65)), linetype = "dotdashed", linewidth = 1., legend = "FFP8 fit")

  do isim = 2, 5
     i = coop_random_index(nsims)
     call fig%curve(r, diffmin(:, i), color = trim(coop_asy_gray_color(0.2)), linetype = "dashed", linewidth = 0.5)
     do j = 0, n
        call coop_chebeval(ncut, 0.d0, rsq(n), vecmin(:, i), rsq(j), fitdiff(j))
     enddo
     call fig%curve(r, fitdiff, color = trim(coop_asy_gray_color(0.65)), linetype = "dotdashed", linewidth = 1.)
     
  enddo
  call coop_asy_legend(fig, "N", 2)
  call fig%close()
  call fig%open(trim(prefix)//"powercut"//COOP_STR_OF(ncut)//"_"//COOP_STR_OF(nsims)//"sims_distr_m"//COOP_STR_OF(m_want)//".txt")
  call fig%init(xlabel = "$c_0 (\mu K)$", ylabel = "$c_1 (\mu K /\mathrm{rad}^2)$")
  call fig%dots(vecmin(1, 1:nsims), vecmin(2, 1:nsims), "black")
  call fig%dot(vecmin(1, 0), vecmin(2, 0), "red", "$\Delta$")
  call fig%expand(0.01, 0.01, 0.01, 0.11)
  call fig%label("$\Delta$   Planck", 0.05, 0.95, "red", alignment="right")
  call fig%label("$\bullet$   FFP8", 0.05, 0.9, color="black", alignment="right")
  call fig%close()
contains

  subroutine fr2vec(d, v)
    COOP_REAL d(0:n)
    COOP_REAL v(ncut)
    COOP_REAL tmpd(n)
    if(m_want .eq. 0)then
       call coop_chebfit(n+1, rsq, d, ncut, 0.d0, rsq(n), v)
    else
       tmpd = d(1:n)/r(1:n)**m_want
       call coop_chebfit(n, rsq(1:n), tmpd, ncut, 0.d0, rsq(n), v)
    endif
  end subroutine fr2vec
  
end program test
