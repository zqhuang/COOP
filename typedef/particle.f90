module coop_particle_mod
  use coop_constants_mod
  use coop_basicutils_mod
  implicit none

#include "constants.h"

  private 

  public:: coop_fermion_get_lnrho, coop_fermion_get_lnp,coop_fermion_get_dlnrho, coop_fermion_get_dlnp, coop_boson_get_lnrho, coop_boson_get_lnp,coop_boson_get_dlnrho, coop_boson_get_dlnp, coop_fermion_get_lnam, coop_boson_get_lnam, coop_fermion_int_q5, coop_fermion_int_q4, coop_fermion_int_q3, coop_fermion_int_kernel5, coop_fermion_int_kernel4, coop_fermion_int_kernel3

  integer,parameter::dl = coop_real_length

  COOP_REAL,dimension(3),parameter::coop_fermion_int_q3 = (/ 0.91320071623_dl, 3.37517484785_dl, 7.7918364724_dl /)
  COOP_REAL,dimension(3),parameter::coop_fermion_int_kernel3 =  (/ 0.0120967050642_dl, 0.583286120914_dl, 0.404617174022_dl /) 

  COOP_REAL,dimension(4),parameter::coop_fermion_int_q4  = (/ 0.4287005867_dl, 2.003741162_dl, 4.783292692_dl,  9.618201254_dl /) 
  COOP_REAL,dimension(4),parameter::coop_fermion_int_kernel4 = (/ 3.720048156e-4_dl, 0.1557779154_dl, 0.6817434900_dl, 0.1621065897_dl /) 

  COOP_REAL,dimension(5),parameter::coop_fermion_int_q5 = (/ 0.3330276475_dl, 1.596709278_dl, 3.756440955_dl, 7.189819223_dl, 12.71011434_dl /) 
  COOP_REAL,dimension(5),parameter::coop_fermion_int_kernel5 =  (/ 1.087668423e-4_dl,  0.06730171005_dl, 0.5294908730_dl, 0.3781574532_dl,  0.02494119690_dl /) 


contains


  subroutine coop_fermion_get_lnam(lnrho, lnam)
    COOP_REAL lnam, lnrho, lower, upper, lr, invam2
    if(lnrho .le. 1.e-2_dl)then
       lnam =  log(lnrho*(1.d0+lnrho/2.d0)/(5_dl/7_dl/coop_pi**2))/2.d0
       return
    endif
    if(lnrho .gt. 1.e2_dl)then
       lnam =  lnrho -  log(coop_Riemannzeta3 * 180_dl/7_dl/coop_pi**4)
       invam2 = exp(-2_dl*lnam)
       lnam = lnrho -  log((coop_Riemannzeta3 * 180_dl/7_dl/coop_pi**4) +( (1350_dl/7_dl*coop_Riemannzeta5/coop_pi**4) + (-6075_dl/4_dl/coop_pi**4 * coop_Riemannzeta7+(172125_dl* coop_Riemannzeta9/8_dl/ coop_pi**4)*invam2)*invam2)*invam2)
       return
    endif
    lower = log((exp(lnrho)-1_dl)*(7_dl*coop_pi**2/5_dl))/2_dl - 1.d-8
    upper = lnrho - log(coop_Riemannzeta3 * 180_dl/7_dl/coop_pi**4) + 1.d-8
    do while(upper - lower .gt. 1.d-7)
       lnam = (lower + upper)/2.d0
       call coop_fermion_get_lnrho(lnam, lr)
       if(lr .gt. lnrho)then
          upper = lnam
       else
          lower = lnam
       endif
    enddo
  end subroutine coop_fermion_get_lnam

!!coop_fermion energy density and pressure, the unit is the energy density of a massless coop_fermion with the same termperature (rho_massless = 7 pi^4/120 T^4)
!!input lnam is ln(mc^2/(kT))
!!output lnrho is ln (rho/rho_massless), ln p is ln(p/rho_massless)
!!in the limit lnam -> - infinity, apparently you get lnrho = 0 and lnp = ln(1/3)
  subroutine coop_fermion_get_lnrho(lnam, lnrho)
    COOP_REAL invam2, corr, lnam, lnrho
    COOP_REAL,dimension(20),parameter::crho = (/ 0.89293935180142925_dl, 1.3851588735560008_dl, 0.59201638940401391_dl, 5.52167768753157942e-2_dl, -6.84994561924834600e-2_dl, -1.78898251764955940e-2_dl, 1.31790803821667645e-2_dl, 5.84725640687740623e-3_dl, -2.65218714293024649e-3_dl, -1.84616659288402659e-3_dl, 4.82404659741137831e-4_dl, 5.58642583176972569e-4_dl, -6.26401384210129567e-5_dl, -1.61520019162570550e-4_dl, -3.30934261743216566e-6_dl, 4.47917623670222794e-5_dl, 6.15434233659787645e-6_dl, -1.15506405054980366e-5_dl, -3.87591649596476188e-6_dl, 3.18654185509784438e-6_dl /)
    if(lnam .le. -2.49_dl)then
       lnrho = log(1_dl + (5_dl/7_dl/coop_pi**2)*exp(2_dl*lnam) )
       return
    endif
    if(lnam .ge. 3.99_dl)then
       invam2 = exp(-2_dl*lnam)
       lnrho = lnam + log((coop_Riemannzeta3 * 180_dl/7_dl/coop_pi**4) +( (1350_dl/7_dl*coop_Riemannzeta5/coop_pi**4) + (-6075_dl/4_dl/coop_pi**4 * coop_Riemannzeta7+(172125_dl* coop_Riemannzeta9/8_dl/ coop_pi**4)*invam2)*invam2)*invam2)
       return
    endif
    call coop_get_cheb_value(20, crho, (lnam - 0.75_dl)/(3.25_dl), lnrho)
  end subroutine coop_fermion_get_lnrho


  subroutine coop_fermion_get_lnp(lnam, lnp)
    COOP_REAL invam2, lnam, lnp, corr
    COOP_REAL,dimension(20),parameter::cp = (/  -1.8856808011378832_dl, -1.2325642746845389_dl, -0.55304532292355002_dl, -8.04535652389012507e-2_dl, 4.86160741219811912e-2_dl, 2.10996354221542545e-2_dl, -5.31380099138448365e-3_dl, -5.02522420384828548e-3_dl, 1.30786399004496982e-4_dl, 1.00515162172776272e-3_dl, 1.70018903096916962e-4_dl, -1.53071349645732708e-4_dl, -6.33387811490677483e-5_dl, 1.22421763460544607e-5_dl, 1.36045554994012240e-5_dl, 1.73976481843332439e-6_dl, -1.66578804580377724e-6_dl, -8.83571732180485515e-7_dl, -1.15092657356077256e-7_dl, 1.54021206255242797e-7_dl /)
    if(lnam .le. -2.49_dl)then
       corr = (5_dl/7_dl/coop_pi**2)*exp(2_dl*lnam)
       lnp = log((1_dl - corr)/3_dl)
       return
    endif
    if(lnam .ge. 3.99_dl)then
       invam2 = exp(-2_dl*lnam)
       lnp = log(((900_dl/7_dl/coop_pi**4*coop_Riemannzeta5) + ((-2025_dl*coop_Riemannzeta7/coop_pi**4) +(172125_dl*coop_RiemannZeta9/2_dl/coop_pi**4)*invam2)*invam2)) - lnam
       return
    endif
    call coop_get_cheb_value(20, cp, (lnam - 0.75_dl)/(3.25_dl), lnp)
  end subroutine coop_fermion_get_lnp

  !!compute d ln rho / d ln am
  subroutine coop_fermion_get_dlnrho(lnam, dlnrho)
    COOP_REAL lnam, dlnrho, corr, invam2
    COOP_REAL,dimension(20),parameter::cdrho = (/ 0.45853294748427892_dl, 0.59814890542859578_dl, 6.46604723315440388e-2_dl, -0.13048676489827551_dl, -3.72780683875114996e-2_dl, 3.81270384863884870e-2_dl, 1.77677707709482345e-2_dl, -1.05345886610559315e-2_dl, -7.42005947098849181e-3_dl, 2.52170510469850589e-3_dl, 2.80538890088544331e-3_dl, -4.47876939099481114e-4_dl, -9.75418076177981598e-4_dl, 1.33083589852769130e-5_dl, 3.17887727370739473e-4_dl, 3.97485986379130440e-5_dl, -9.38232585120578154e-5_dl, -2.40406623904638736e-5_dl, 2.99526741850547889e-5_dl, 1.34980095732365434e-5_dl /)
    if(lnam .le. -2.49_dl)then
       corr = (5_dl/7_dl/coop_pi**2)*exp(2_dl*lnam)
       dlnrho = 2_dl*corr/(1_dl + corr)
       return
    endif
    if(lnam .gt. 3.99_dl)then
       invam2 = exp(-2_dl*lnam)
       dlnrho = 1_dl + ( (-15_dl*coop_Riemannzeta5/coop_Riemannzeta3) &
            + (45_dl*(10_dl*coop_Riemannzeta5**2 + 21_dl*coop_Riemannzeta3*coop_Riemannzeta7)/4_dl/coop_Riemannzeta3**2)*invam2)*invam2
       return
    endif
    call coop_get_cheb_value(20, cdrho, (lnam-0.75_dl)/3.25_dl, dlnrho)
  end subroutine coop_fermion_get_dlnrho



    !!compute d ln rho / d ln am
  subroutine coop_fermion_get_dlnp(lnam, dlnp)
    COOP_REAL lnam, dlnp, corr, invam2
    COOP_REAL, dimension(20),parameter::cdrho = (/ -0.42955931767638317_dl, -0.57929706031950723_dl, -0.10061753508327952_dl, 0.10137410858790824_dl, 0.47912127913089325e-1_dl , -0.18296247083641887e-1_dl , -0.17009826170227303e-1_dl , 0.13239071953237954e-2_dl , 0.46372977504751667e-2_dl , 0.67998892886111039e-3_dl , -0.92969391416750698e-3_dl , -0.36634500130460678e-3_dl , 0.10647905926389047e-3_dl , 0.10129351444516322e-3_dl , 0.85424263264848575e-5_dl , -0.16063650866938333e-4_dl , -0.75048662407985448e-5_dl , 0.98439667538101790e-7_dl, 0.17632977601903781e-5_dl , 0.97383174483000743e-6_dl  /)
    if(lnam .le. -2.49_dl)then
       corr = (5_dl/7_dl/coop_pi**2)*exp(2_dl*lnam)
       dlnp = -2.* corr/(1.-corr)
       return
    endif
    if(lnam .gt. 3.99_dl)then
       invam2 = exp(-2_dl*lnam)
       dlnp = - ((-2.*2025_dl*coop_Riemannzeta7/coop_pi**4) +(4.*172125_dl*coop_RiemannZeta9/2_dl/coop_pi**4)*invam2) *invam2 &
            / (((900_dl/7_dl/coop_pi**4*coop_Riemannzeta5) + ((-2025_dl*coop_Riemannzeta7/coop_pi**4) +(172125_dl*coop_RiemannZeta9/2_dl/coop_pi**4)*invam2)*invam2)) &
            - 1._dl
       return
    endif
    call coop_get_cheb_value(20,cdrho, (lnam-0.75_dl)/3.25_dl, dlnp)
  end subroutine coop_fermion_get_dlnp
!!------------------------ Coop_boson -----------------------------------------

  subroutine coop_boson_get_lnam(lnrho, lnam)
    COOP_REAL lnam, lnrho, lower, upper, lr, invam2
    if(lnrho .le. 1.d-2)then
       lnam = log (lnrho * (1.d0 + lnrho/2.d0) /  (5_dl/4_dl/coop_pi**2))/2.d0
       return
    endif
    if(lnrho .ge. 1.d2)then
       lnam = lnrho - log((coop_Riemannzeta3 * 30_dl/coop_pi**4))
       invam2 = exp(-2_dl*lnam)
       lnam = lnrho - log((coop_Riemannzeta3 * 30_dl/coop_pi**4) +( (180_dl*coop_Riemannzeta5/coop_pi**4) + (-1350_dl/coop_pi**4 * coop_Riemannzeta7 + (37800_dl * coop_Riemannzeta9/ coop_pi**4)*invam2)*invam2)*invam2)
       return
    endif
    lower = log((exp(lnrho)-1_dl)*(4_dl*coop_pi**2/5_dl))/2_dl - 1.d-8
    upper = lnrho - log(coop_Riemannzeta3 * 30_dl/coop_pi**4) + 1.d-8
    do while(upper - lower .gt. 1.d-7)
       lnam = (lower + upper)/2.d0
       call coop_boson_get_lnrho(lnam, lr)
       if(lr .gt. lnrho)then
          upper = lnam
       else
          lower = lnam
       endif
    enddo
  end subroutine coop_boson_get_lnam




!!coop_boson energy density and pressure, the unit is the energy density of a massless bonson with the same termperature rho_massless = pi^4/15 T^4
!!input lnam is ln(mc^2/(kT))
!!output lnrho is ln (rho/rho_massless), ln p is ln(p/rho_massless)
!!in the limit lnam -> - infinity, apparently you get lnrho = 0 and lnp = ln(1/3)
  subroutine coop_boson_get_lnrho(lnam, lnrho)
    COOP_REAL invam2, corr, lnam, lnrho
    COOP_REAL,dimension(20),parameter::crho = (/ 0.85920499059112343_dl, 1.3962721264359432_dl, 0.70432294674814611_dl, 0.14297095706190513_dl, -6.55811511010368081e-2_dl, -4.63591428381798787e-2_dl, 5.01875599077331383e-3_dl, 1.41690006104390458e-2_dl, 2.21177480941542586e-3_dl, -3.70703629807588638e-3_dl, -1.66893509068491972e-3_dl, 7.15614414586352820e-4_dl, 7.23143566964830971e-4_dl, -3.35664246204683407e-5_dl, -2.42009019038368436e-4_dl, -5.66482185018125605e-5_dl, 6.35841417332966958e-5_dl, 3.44090360869310461e-5_dl, -1.09428321341706716e-5_dl, -1.68472313639487885e-5_dl /)
    if(lnam .le. -3.99_dl)then
       corr = (5_dl/4_dl/coop_pi**2)*exp(2_dl*lnam)
       lnrho = log(1_dl + corr)
       return
    endif
    if(lnam .ge. 3.99_dl)then
       invam2 = exp(-2_dl*lnam)
       lnrho = lnam + log((coop_Riemannzeta3 * 30_dl/coop_pi**4) +( (180_dl*coop_Riemannzeta5/coop_pi**4) + (-1350_dl/coop_pi**4 * coop_Riemannzeta7 + (37800_dl * coop_Riemannzeta9/ coop_pi**4)*invam2)*invam2)*invam2)
       return
    endif
    call coop_get_cheb_value(20, crho,  lnam/4._dl, lnrho)
  end subroutine coop_boson_get_lnrho

  subroutine coop_boson_get_lnp(lnam, lnp)
    COOP_REAL invam2, corr, lnam, lnp
    COOP_REAL,dimension(20),parameter::cp = (/ -1.8285425606893564_dl, -1.1995994663208203_dl, -0.63579323906066232_dl, -0.16456435145041079_dl, 3.22461274266371070e-2_dl, 4.02357879959531692e-2_dl, 5.56128174481804183e-3_dl, -7.42808998219154047e-3_dl, -3.68630245376660396e-3_dl, 5.37511637030240472e-4_dl, 1.03818428724796819e-3_dl, 2.31230717643599081e-4_dl, -1.53291700135450192e-4_dl, -1.07074751727548019e-4_dl, -6.96618941049786799e-6_dl, 2.02895088393592875e-5_dl, 1.00029976391815033e-5_dl, -7.29511466983895998e-8_dl, -2.23609786187199272e-6_dl, -1.23155707342790243e-6_dl /)
    if(lnam .le. -3.99_dl)then
       corr = (5_dl/4_dl/coop_pi**2)*exp(2_dl*lnam)
       lnp = log((1_dl - corr)/3_dl)
       return
    endif
    if(lnam .ge. 3.99_dl)then
       invam2 = exp(-2_dl*lnam)
       lnp = log(((120_dl/coop_pi**4*coop_Riemannzeta5) + ((-1800_dl*coop_Riemannzeta7/coop_pi**4) +(75600_dl*coop_RiemannZeta9/coop_pi**4)*invam2)*invam2)) - lnam       
       return
    endif
    call coop_get_cheb_value(20, cp, lnam/4._dl, lnp)
  end subroutine coop_boson_get_lnp
  

  !!compute d ln rho / d ln am
  subroutine coop_boson_get_dlnrho(lnam, dlnrho)
    COOP_REAL lnam, dlnrho, corr, invam2
    COOP_REAL,dimension(20),parameter::cdrho = (/ 0.41655827483174224_dl, 0.59177723567654272_dl, 0.13498019867091371_dl, -0.11254573750954061_dl, -7.94769037761112740e-2_dl, 1.86165504389742276e-2_dl, 3.64196607223354896e-2_dl, 3.56021952692404855e-3_dl, -1.31738437483698515e-2_dl, -5.28694710409113346e-3_dl, 3.50481795939059987e-3_dl, 3.05761662838206988e-3_dl, -4.35482712001917550e-4_dl, -1.28139839893068341e-3_dl, -2.23845853250542184e-4_dl, 4.12439057719122951e-4_dl, 1.91001532358664038e-4_dl, -9.65718422490377602e-5_dl, -1.18271863089807312e-4_dl, 1.35735494059749740e-6_dl /)
    if(lnam .le. -2.49_dl)then
       corr = (5_dl/4_dl/coop_pi**2)*exp(2_dl*lnam)
       dlnrho = 2_dl*corr/(1_dl + corr)
       return
    endif
    if(lnam .gt. 3.99_dl)then
       invam2 = exp(-2_dl*lnam)
       dlnrho = 1_dl + ( (-12_dl*coop_Riemannzeta5/coop_Riemannzeta3) &
            + ((72_dl*coop_Riemannzeta5**2 + 180_dl*coop_Riemannzeta3*coop_Riemannzeta7)/coop_Riemannzeta3**2)*invam2)*invam2
       return
    endif
    call coop_get_cheb_value(20, cdrho, lnam/4._dl, dlnrho)
  end subroutine coop_boson_get_dlnrho

  subroutine coop_boson_get_dlnp(lnam, dlnp)
    COOP_REAL lnam, dlnp, corr, invam2
    COOP_REAL, dimension(20),parameter::cdrho =   (/ -0.38445753351392464_dl, -0.56507922513479980_dl, -0.16911534555220517_dl, 0.70714005528823232e-1_dl , 0.77731140160035683e-1_dl , 0.62217410104989388e-2_dl , -0.22858414171934623e-1_dl , -0.10462112478909665e-1_dl , 0.31397644217964569e-2_dl , 0.42830853648848930e-2_dl , 0.72075785662759496e-3_dl , -0.90785063326038519e-3_dl , -0.55131656392327173e-3_dl , 0.11882642415894909e-4_dl , 0.14421947573201171e-3_dl , 0.60608967122437350e-4_dl , -0.86129681024491393e-5_dl , -0.19507626000606408e-4_dl , -0.90584422331560452e-5_dl , 0.43136694416906456e-6_dl  /)
    if(lnam .le. -3.99_dl)then
       corr = (5_dl/4_dl/coop_pi**2)*exp(2_dl*lnam)
       dlnp = -2.* corr/(1.-corr)
       return
    endif
    if(lnam .gt. 3.99_dl)then
       invam2 = exp(-2_dl*lnam)
       dlnp = -  (((-2.*1800_dl*coop_Riemannzeta7/coop_pi**4) +(4.*75600_dl*coop_RiemannZeta9/coop_pi**4)*invam2)*invam2)  &
            / ((120_dl/coop_pi**4*coop_Riemannzeta5) + ((-1800_dl*coop_Riemannzeta7/coop_pi**4) +(75600_dl*coop_RiemannZeta9/coop_pi**4)*invam2)*invam2) &
            - 1._dl
       return
    endif
    call coop_get_cheb_value(20, cdrho, lnam/4._dl, dlnp)
  end subroutine coop_boson_get_dlnp

end module coop_particle_mod


