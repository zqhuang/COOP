program bgtest
  use coop_wrapper_background
  implicit none
#include "constants.h"  
  type(coop_cosmology_background)::bg, bg_save
  COOP_REAL:: Qcpl, tracking_n, dlnQdphi, dUdphi, d2Udphi2, omegab, omegac, h, omegar, omegam
  COOP_INT:: iQ, in, index_DE, index_CDM, i
  COOP_INT,parameter::nn = 20, nQ = 15
  COOP_REAL,parameter::apiv = 0.1d0
  COOP_REAL::deform(0:nn, 0:nQ), nlist(0:nn)
  type(coop_asy)::fig
  dlnQdphi = 0.d0
  dUdphi = 0.d0
  d2Udphi2 = 0.d0
  omegab = 0.049d0
  omegac = 0.24d0
  omegam = omegab + omegac
  call bg_save%init(h=0.68d0)
  call bg_save%add_species(coop_baryon(omegab))
  call bg_save%add_species(coop_radiation(bg%Omega_radiation()))
  call bg_save%add_species(coop_neutrinos_massless(bg%Omega_massless_neutrinos_per_species()*(bg%Nnu())))
  index_CDM = 4
  index_DE = 5
  omegar = bg_save%species(2)%omega+bg_save%species(3)%omega
  call fig%open("deform.txt")
  call fig%init(xlabel="$n$", ylabel="$\log_{10}Q$")
  do iQ = 0, nQ
     Qcpl = 10.d0**((iQ-nQ)*2.d0/nQ)
     write(*,*) Qcpl
     do in = 0, nn
        tracking_n = 0.05 + 2.d0/nn * in
        nlist(in) = tracking_n
        bg = bg_save
        call coop_background_add_coupled_DE(bg, Omega_c = omegac, Q = Qcpl, tracking_n = tracking_n, dlnQdphi = dlnQdphi, dUdphi = dUdphi, d2Udphi2 = d2Udphi2)
        deform(in, iQ) = ((1.d0 + bg%species(index_DE)%pa2(apiv)/( bg%species(index_DE)%rhoa2(apiv) + bg%species(index_CDM)%rhoa2(apiv)-(3.d0*omegac)/apiv)) - (tracking_n/(tracking_n+2.d0)))
     enddo
  enddo
  call fig%density(deform, 0.05d0, 2.05d0, -3.d0, 0.d0, "$\epsilon/\epsilon_0$")
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
