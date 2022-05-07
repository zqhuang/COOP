program test
  use coop_wrapper_utils
  use coop_background_mod
  implicit none
#include "constants.h"
  type(coop_cosmology_background) bg
  integer,parameter::n=740, nw=2, nbins = 7
  type(coop_file) fp
  type(coop_asy) fig
  real, parameter::ymax = 0.12, ymin = -0.1
  COOP_REAL,parameter::step = 0.05
  COOP_LONG_STRING line
  COOP_SHORT_STRING name
  COOP_REAL w
  COOP_REAL a(n), zcmb(n), zhel(n), mb(n), dmb(n), dz, mu_theory(n, -nw:nw), Mabs, junk
  COOP_REAL, parameter::alpha = 0.141, beta = 3.10, intr_dm = 0.13
  integer i, iw, ind(n), set(n), nlabel
  COOP_REAL binned_z(nbins), binned_mu(nbins), binned_dmu(nbins), weight(nbins), strech, color, dstrech, dcolor, cov_ms, cov_mc, cov_sc, thirdvar, dthirdvar
  
  call fp%open("../data/jla_lcparams.txt","r")
  read(fp%unit, "(A)") line
  do i=1, n
     read(fp%unit, *) name, zcmb(i), zhel(i), dz, mb(i), dmb(i), strech, dstrech, color, dcolor, thirdvar, dthirdvar, cov_ms, cov_mc, cov_sc, set(i)
     mb(i) = mb(i) - 5.d0*log10((1.+zhel(i))/(1.+zcmb(i))) + alpha*(strech-1.d0) - beta*color
     dmb(i) = sqrt(dmb(i)**2 + (5.d0*log10(1.d0+0.001/zcmb(i)))**2 + intr_dm**2 + alpha**2*dstrech**2 + beta**2*dcolor**2 + 2.d0*alpha*cov_ms - 2.d0*beta*cov_mc - 2.d0*alpha*beta*cov_sc)
  enddo
  call fp%close()
  call coop_quicksortacc(zcmb, ind)
  mb = mb(ind)
  dmb = dmb(ind)
  a = 1.d0/(1.d0+zcmb)
  do iw=-nw, nw
     w =-1+ step*iw
     call bg%init(h=COOP_REAL_OF(0.682d0))
     call bg%add_species(coop_baryon(COOP_REAL_OF(0.047d0)))
     call bg%add_species(coop_cdm(COOP_REAL_OF(0.253d0)))
     call bg%add_species(coop_radiation(bg%Omega_radiation()))
     call bg%add_species(coop_neutrinos_massive( bg%Omega_nu_from_mnu_eV(0.06d0),bg%Omega_massless_neutrinos_per_species()))
     call bg%add_species(coop_neutrinos_massless(bg%Omega_massless_neutrinos_per_species()*2.))
     call bg%add_species(coop_de_w0(bg%Omega_k(), w))
   !  call bg%add_species(coop_de_quintessence(bg%Omega_k(), epsilon_s, 0.d0, 0.d0))
     call bg%setup_background()
     do i=1, n
        mu_theory(i, iw) = 5.d0*log10(bg%luminosity_distance(a(i))/bg%H0Mpc())+25.
     enddo
     call bg%free()
  enddo
  Mabs = sum((mu_theory(:,0)-mb)/dmb**2)/sum(1.d0/dmb**2)
  print*, "M = ", Mabs

  mb = mb + Mabs - mu_theory(:, 0)
  binned_z = 0.d0
  binned_mu = 0.d0
  weight = 0.d0
  binned_dmu = 0.d0
  do i=1, n
     iw = ceiling((zcmb(i) + 0.1d0)/0.2d0)
     if(iw .lt. 1 .or. iw.gt. nbins)stop "bin range overflow"
     weight(iw) = weight(iw) + 1.d0/dmb(i)**2
     binned_z(iw) = binned_z(iw) + zcmb(i)/dmb(i)**2
     binned_mu(iw) = binned_mu(iw) + mb(i)/dmb(i)**2
     binned_dmu(iw) = binned_dmu(iw) + 1.d0/dmb(i)**2
  enddo
  if(any(weight .le. 0.d0)) stop "empty bin"
  binned_mu = binned_mu/weight
  binned_z = binned_z/weight
  binned_dmu = 1.d0/sqrt(binned_dmu)
  do iw=-nw, nw
     if(iw.ne.0) &
          mu_theory(:, iw) = mu_theory(:, iw) - mu_theory(:, 0)
  enddo
  mu_theory(:, 0) = 0.
  call fig%open("JLAsupernova_w.txt")
  call fig%init(xlabel = "$z$", ylabel = "$\Delta\mu$", xmin = 0., xmax = 1.35, ymin = ymin, ymax = ymax, doclip = .true.)
  nlabel = n-38
  do iw = -nw, nw
     call coop_asy_interpolate_curve(fig, xraw = zcmb, yraw = mu_theory(:, iw), interpolate = "LinearLinear", color=fig%color(iw+nw+1), linetype = fig%linetype(iw+nw+1))
     call coop_asy_label(fig, x = real(zcmb(nlabel)), y = real(mu_theory(nlabel, iw))+0.005, label = "$w="//trim(coop_num2str(-1.d0+step*iw))//"$")
  enddo
  do i=1, nbins
     call coop_asy_error_bar(fig, x = binned_z(i), y = binned_mu(i), dy_plus = binned_dmu(i), dy_minus = binned_dmu(i))
  enddo
  call coop_asy_label(fig, x = 0.22, y = ymin+0.03, label = "$\Omega_m = 0.3$")
  call coop_asy_label(fig, x = 0.5, y = ymax-0.008, label = "JLA SN data, $\Lambda$CDM model subtracted")
  call fig%close()
end program test
