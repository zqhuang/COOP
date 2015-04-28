program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  COOP_REAL::tracking_n = 0.3d0
  COOP_REAL::h0 = 0.68d0
  COOP_REAL::omega_b = 0.049d0
  COOP_REAL::omega_c = 0.265d0
  type(coop_cosmology_firstorder)::fod

  COOP_INT,parameter::nk = 4, nz = 256
  COOP_REAL,dimension(nk)::kMpc = (/ 0.001, 0.005, 0.01, 0.1 /)
  COOP_REAL,dimension(nk)::k
  COOP_REAL,dimension(nz)::z, growth
  type(coop_asy)::fig
  COOP_INT:: ik, iz, iQ
  COOP_REAL::Q
  call coop_set_uniform(nz, z, 0.d0, 2.5d0)
  do iQ = -1, 3
     Q  = iQ * 0.1d0
     if(Q.ge.0.d0)then
        call fig%open("growth_Q0p"//COOP_STR_OF(iQ)//".txt")        
     else
        call fig%open("growth_LCDM.txt")
     endif
     call fig%init(xlabel = "$z$", ylabel = "$\Phi(z; k)$",  caption="Coupled Dark Energy $Q = \frac{\partial \ln m_{DM}}{\partial \phi}$", ymin = 0.75, ymax = 1.05)

     if(Q.ge. 0.d0)then
        call fod%set_standard_cosmology(Omega_b=omega_b, Omega_c= omega_c, h = h0, tau_re = 0.08d0, As = 2.2d-9, ns = 0.96d0, de_Q = Q, de_tracking_n = tracking_n)
     else
        call fod%set_standard_cosmology(Omega_b=omega_b, Omega_c= omega_c, h = h0, tau_re = 0.08d0, As = 2.2d-9, ns = 0.96d0)
     endif
     k = kMpc/fod%H0Mpc() !!convert k to H0 unit
     call fod%compute_source(0)

     do ik= 1, nk
        !$omp parallel do
        do iz = 1, nz
           growth(iz) = fod%growth_of_z(z = z(iz), k = k(ik))
        enddo
        !$omp end parallel do
        call fig%curve(z, growth, linewidth = fig%linewidth(ik), color = fig%color(ik), linetype = fig%linetype(ik), legend="$k = "//COOP_STR_OF(kMpc(ik))//" {\rm Mpc}^{-1}$")     
     enddo
     if(Q.ge. 0.d0)then
        call fig%label("$V = V_0/\phi^{"//COOP_STR_OF(tracking_n)//"},\, Q="//COOP_STR_OF(Q)//"$", 0.45, 0.6)
     else
        call fig%label("$\Lambda$CDM", 0.45, 0.6)        
     endif
     call fig%legend(0.45, 0.35)


     call fig%close()
  enddo
end program test
