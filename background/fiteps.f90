program bgtest
  use coop_wrapper_background
  implicit none
#include "constants.h"  
  type(coop_cosmology_background)::bg, bg_save
  COOP_INT, parameter::n = 256
  COOP_REAL:: Qcpl, tracking_n, dlnQdphi, dUdphi, d2Udphi2, omegab, omegac, h, omegar, omegam
  COOP_REAL:: a(n), eps(n), omde(n), omm(n), omr(n), z(n), epss(n), epsfit(n), rhoa3, aeq, epseq, epsilon_infty, zeta_s, epsilon_s, lambda
  COOP_INT::i, iU, j, iearly, index_DE, index_CDM
  type(coop_asy)::fig
  tracking_n = 1.d0 
  dlnQdphi = 0.d0
  dUdphi = 0.d0
  d2Udphi2 = 0.d0
  omegab = 0.049d0
  omegac = 0.265d0
  omegam = omegab + omegac
  call bg_save%init(h=0.68d0)
  call bg_save%add_species(coop_baryon(omegab))
  call bg_save%add_species(coop_radiation(bg%Omega_radiation()))
  call bg_save%add_species(coop_neutrinos_massless(bg%Omega_massless_neutrinos_per_species()*(bg%Nnu())))
  index_CDM = 4
  index_DE = 5
  omegar = bg_save%species(2)%omega+bg_save%species(3)%omega

  call coop_set_uniform(n, a, 1.d-2, 1.d0, logscale = .true.)
  z = 1.d0/a-1.d0

  call fig%open("epsa.txt")
  call fig%init(xlabel="$a$", ylabel="$\epsilon$", width = 8., height = 6., xmin = real(a(1)), xmax =real(a(n)), ymin = -0.01, ymax = 0.2, doclip = .true.)
  Qcpl = 0.d0
  do iU = 0, 0
     dUdphi = -iU*0.2d0 ! + tracking_n 
     bg = bg_save
     call coop_background_add_coupled_DE(bg, Omega_c = omegac, Q = Qcpl, tracking_n = tracking_n, dlnQdphi = dlnQdphi, dUdphi = dUdphi, d2Udphi2 = d2Udphi2)
     aeq = 0.d0
     iearly = 0
     do i = 1, n
        rhoa3 = bg%rhoa4(a(i))/a(i)/3.d0
        omm(i) = omegam/rhoa3
        omr(i) = omegar/a(i)/rhoa3
        omde(i) = 1.d0 - omm(i) - omr(i)
        eps(i) = -bg%HdotbyHsq(a(i))
        epss(i) = eps(i) - 1.5d0*omm(i) - 2.d0*omr(i)
        write(*,"(3E14.5)") a(i), 1.d0+bg%species(index_DE)%pa2(a(i))/(bg%species(index_DE)%rhoa2(a(i)) + bg%species(index_CDM)%rhoa2(a(i))- 3.d0*omegac/a(i)),  1.d0+bg%species(index_DE)%pa2(a(i))/(bg%species(index_DE)%rhoa2(a(i)))
        if(omde(i) .gt. 0.5d0 .and. aeq.eq. 0.d0 .and. i.gt.1)then
           aeq = (a(i-1)*(omde(i)-0.5d0)+a(i)*(0.5d0-omde(i-1)))/(omde(i)-omde(i-1))

           epseq = (epss(i-1)*(omde(i)-0.5d0)+epss(i)*(0.5d0-omde(i-1)))/(omde(i)-omde(i-1))
        endif
        if(iearly.eq.0 .and. omde(i) .gt. 0.05d0)then
           iearly = i
        endif
     enddo
     stop
     epsilon_infty = epss(iearly)/omde(iearly)
     epsilon_s = ((sqrt(epss(n)/omde(n)) - sqrt(epsilon_infty)) / (quint_f(1.d0/aeq)) + sqrt(2.d0*epsilon_infty))**2
     zeta_s = max(min( ((sqrt(epseq*2.d0) - sqrt(epsilon_infty))/(sqrt(epsilon_s) - sqrt(2.d0*epsilon_infty))-quint_f(1.d0))/quint_f2(1.d0), 1.d0), -1.d0)
     print*, a(iearly), aeq, epsilon_s, epsilon_infty, zeta_s
     do i=1, n
        epsfit(i) = omde(i)*(sqrt(epsilon_infty) + (sqrt(epsilon_s)-sqrt(2.d0*epsilon_infty))*(quint_f(a(i)/aeq)+zeta_s*quint_f2(a(i)/aeq)))**2
     enddo
     call coop_asy_curve(fig, x = a, y = epss, linetype = fig%linetype(iU+1), color=fig%color(iU+1), legend = "$Q="//COOP_STR_OF(Qcpl)//", dU/d\phi="//COOP_STR_OF(dUdphi)//"$", linewidth = 1.8)
     call coop_asy_curve(fig, x = a, y = epsfit, linetype = fig%linetype(iU+1), color=fig%color(iU+1))
  enddo
  call fig%legend(0.1, 0.9)
  call fig%close()

contains

  function quint_f(x)
    COOP_REAL::x, quint_f
    quint_f = sqrt(1.d0/x**3+1.d0)-log(x**1.5d0+sqrt(1.d0+x**3))/x**3
  end function quint_f

  function quint_f2(x)
    COOP_REAL::x, quint_f2
    quint_f2 = coop_sqrt2*(1.d0-log(1.d0+x**3)/x**3) - quint_f(x)
  end function quint_f2

    
  
end program bgtest
