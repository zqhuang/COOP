module coop_halofit_mod
  use coop_wrapper_firstorder
  implicit none

#include "constants.h"

  private


  public::coop_halofit_get_power

contains

  !!return the dimensionless P
  subroutine coop_halofit_get_power(cosmology, z, nk, k, Pl, Pnl)
    class(coop_cosmology_firstorder)::cosmology
    COOP_INT nk
    COOP_REAL k(nk), Pnl(nk), Pl(nk), z
    COOP_REAL om_m, om_v, fnu, w0
    integer itf
    COOP_REAL a,rk
    COOP_REAL sig,rknl,rneff,rncur,d1,d2
    COOP_REAL diff, rup, rdown, rmid
    integer i, nloops
    call cosmology%get_matter_power(z, nk, k, pl)
    a = 1.d0/(1.d0+z)
    om_m = cosmology%omega_m*a/cosmology%H2a4(a)
    om_v = 1.d0 - om_m
    fnu = cosmology%Omega_massivenu / Cosmology%Omega_m
    w0 = O0_DE(cosmology)%wofa(coop_scale_factor_today)

    rdown = 0.1d0/k(nk)  
    rup = 1.d0 !!  k ~ H_0
    if(cosmology%sigma_Gaussian_R_quick(z, rdown) .le. 0.97d0)then
       Pnl = Pl
       return
    elseif(cosmology%sigma_Gaussian_R_quick(z, rup) .ge. 1.03d0)then
       call coop_return_error("halofit_power", "linear power overflow", "stop")
    endif
    nloops = 1
    do 
       rmid=sqrt(rdown*rup)
       sig = cosmology%sigma_Gaussian_R_quick(z, rmid)
       if(sig .ge. 1.03d0)then
          rdown = rmid
       elseif(sig .le. 0.97d0)then
          rup = rmid
       else
          exit
       endif
       nloops = nloops+1
       if(nloops .gt. 200) call coop_return_error("halofit_power", "r_G not converging", "stop")
    end do
    nloops = 1
    do 
       rmid=sqrt(rdown*rup)
       sig = cosmology%sigma_Gaussian_R(z, rmid)
       if(sig .ge. 1.003d0)then
          rdown = rmid
       elseif(sig .le. 0.997d0)then
          rup = rmid
       else
          exit
       endif
       nloops = nloops+1
       if(nloops .gt. 200) call coop_return_error("halofit_power", "r_G not converging", "stop")
    end do
    call cosmology%sigma_Gaussian_R_with_dervs(z, rmid, sig, d1, d2)
    rmid = rmid *exp( log(1.d0/sig**2)/d1)
    call cosmology%sigma_Gaussian_R_with_dervs(z, rmid, sig, d1, d2)
    rknl = (1.d5/coop_SI_c)/rmid
    rneff  = -3.d0 - d1
    rncur = - d2
    !$omp parallel do
    do i=1, nk
       if(k(i)*rmid .lt. 1.d-2)then
          Pnl(i) = pl(i)
       else
          call coop_halofit_one(w0, om_m, om_v, fnu, cosmology%omega_m,  k(i)*(1.d5/coop_SI_c) ,rneff,rncur,rknl,pl(i),pnl(i))   ! halo fitting formula 
       endif
    enddo
    !$omp end parallel do
  end subroutine coop_halofit_get_power



  subroutine coop_halofit_one(w0, om_m,om_v,fnu,omm0, rk,rn,rncur,rknl,plin,pnl)
    ! halo model nonlinear fitting formula as described in 
    ! Appendix C of Smith et al. (2002)
    COOP_REAL om_m,om_v,fnu,omm0, w0
    COOP_REAL gam,a,b,c,xmu,xnu,alpha,beta,f1,f2,f3
    COOP_REAL rk,rn,plin,pnl,pq,ph,plinaa
    COOP_REAL rknl,y,rncur
    COOP_REAL f1a,f2a,f3a,f1b,f2b,f3b,frac
    !RT12 Oct: the halofit in Smith+ 2003 predicts a smaller power
    !than latest N-body simulations at small scales.
    !Update the following fitting parameters of gam,a,b,c,xmu,xnu,
    !alpha & beta from the simulations in Takahashi+ 2012.
    !The improved halofit accurately provide the power spectra for WMAP
    !cosmological models with constant w.

    gam=0.1971-0.0843*rn+0.8460*rncur
    a=1.5222+2.8553*rn+2.3706*rn*rn+0.9903*rn*rn*rn+ &
         0.2250*rn*rn*rn*rn-0.6038*rncur+0.1749*om_v*(1.+w0)
    a=10**a      
    b=10**(-0.5642+0.5864*rn+0.5716*rn*rn-1.5474*rncur+ &
         0.2279*om_v*(1.+w0))
    c=10**(0.3698+2.0404*rn+0.8161*rn*rn+0.5869*rncur)
    xmu=0.
    xnu=10**(5.2105+3.6902*rn)
    alpha=abs(6.0835+1.3373*rn-0.1959*rn*rn-5.5274*rncur)
    beta=2.0379-0.7354*rn+0.3157*rn**2+1.2490*rn**3+ &
         0.3980*rn**4-0.1682*rncur + fnu*(1.081 + 0.395*rn**2)

    if(abs(1-om_m).gt.0.01) then ! omega evolution 
       f1a=om_m**(-0.0732)
       f2a=om_m**(-0.1423)
       f3a=om_m**(0.0725)
       f1b=om_m**(-0.0307)
       f2b=om_m**(-0.0585)
       f3b=om_m**(0.0743)       
       frac=om_v/(1.-om_m) 
       f1=frac*f1b + (1-frac)*f1a
       f2=frac*f2b + (1-frac)*f2a
       f3=frac*f3b + (1-frac)*f3a
    else         
       f1=1.0
       f2=1.
       f3=1.
    endif
    y=(rk/rknl)
    ph=a*y**(f1*3)/(1+b*y**(f2)+(f3*c*y)**(3-gam))
    ph=ph/(1+xmu*y**(-1)+xnu*y**(-2))*(1+fnu*(0.977-18.015*(omm0-0.3)))
    plinaa=plin*(1+fnu*47.48*rk**2/(1+1.5*rk**2))
    pq=plin*(1+plinaa)**beta/(1+plinaa*alpha)*exp(-y/4.0-y**2/8.0)
    pnl=pq+ph
  end subroutine coop_halofit_one


end module coop_halofit_mod
