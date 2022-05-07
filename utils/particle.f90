module particle_utils
  use general_utils
  implicit none

  real(dl),dimension(3),parameter::fermion_int_q3 = (/ 0.91320071623d0, 3.37517484785d0, 7.7918364724d0 /)
  real(dl),dimension(3),parameter::fermion_int_kernel3 =  (/ 0.0120967050642_dl, 0.583286120914_dl, 0.404617174022_dl /) 

  real(dl),dimension(4),parameter::fermion_int_q4  = (/ 0.4287005867d0, 2.003741162d0, 4.783292692d0,  9.618201254d0 /) 
  real(dl),dimension(4),parameter::fermion_int_kernel4 = (/ 3.720048156e-4_dl, 0.1557779154_dl, 0.6817434900_dl, 0.1621065897_dl /) 

  real(dl),dimension(5),parameter::fermion_int_q5 = (/ 0.3330276475_dl, 1.596709278_dl, 3.756440955_dl, 7.189819223_dl, 12.71011434_dl /) 
  real(dl),dimension(5),parameter::fermion_int_kernel5 =  (/ 1.087668423e-4_dl,  0.06730171005_dl, 0.5294908730_dl, 0.3781574532_dl,  0.02494119690_dl /) 


contains

!!fermion energy density and pressure, the unit is the energy density of a massless fermion with the same termperature (rho_massless = 7 pi^4/120 T^4)
!!input lnam is ln(mc^2/(kT))
!!output lnrho is ln (rho/rho_massless), ln p is ln(p/rho_massless)
!!in the limit lnam -> - infinity, apparently you get lnrho = 0 and lnp = ln(1/3)
  subroutine fermion_get_lnrho(lnam, lnrho)
    real(dl) invam2, corr, lnam, lnrho
    real(dl),dimension(20),parameter::crho = (/ 0.89293935180142925_dl, 1.3851588735560008_dl, 0.59201638940401391_dl, 5.52167768753157942e-2_dl, -6.84994561924834600e-2_dl, -1.78898251764955940e-2_dl, 1.31790803821667645e-2_dl, 5.84725640687740623e-3_dl, -2.65218714293024649e-3_dl, -1.84616659288402659e-3_dl, 4.82404659741137831e-4_dl, 5.58642583176972569e-4_dl, -6.26401384210129567e-5_dl, -1.61520019162570550e-4_dl, -3.30934261743216566e-6_dl, 4.47917623670222794e-5_dl, 6.15434233659787645e-6_dl, -1.15506405054980366e-5_dl, -3.87591649596476188e-6_dl, 3.18654185509784438e-6_dl /)
    if(lnam .le. -2.49d0)then
       corr = (5.d0/7.d0/const_pi2)*exp(2.d0*lnam)
       lnrho = log(1.d0 + corr)
       return
    endif
    if(lnam .ge. 3.99d0)then
       invam2 = exp(-2.d0*lnam)
       lnrho = lnam + log((const_Riemannzeta3 * 180.d0/7.d0/const_pi4) +( (1350.d0/7.d0*const_Riemannzeta5/const_pi4) + (-6075.d0/4.d0/const_pi4 * const_Riemannzeta7+(172125.d0* const_Riemannzeta9/8.d0/ const_pi4)*invam2)*invam2)*invam2)
       return
    endif
    call get_cheb_value(crho, 20,  (lnam - 0.75d0)/(3.25d0), lnrho)
  end subroutine fermion_get_lnrho


  subroutine fermion_get_lnp(lnam, lnp)
    real(dl) invam2, lnam, lnp, corr
    real(dl),dimension(20),parameter::cp = (/  -1.8856808011378832_dl, -1.2325642746845389_dl, -0.55304532292355002_dl, -8.04535652389012507e-2_dl, 4.86160741219811912e-2_dl, 2.10996354221542545e-2_dl, -5.31380099138448365e-3_dl, -5.02522420384828548e-3_dl, 1.30786399004496982e-4_dl, 1.00515162172776272e-3_dl, 1.70018903096916962e-4_dl, -1.53071349645732708e-4_dl, -6.33387811490677483e-5_dl, 1.22421763460544607e-5_dl, 1.36045554994012240e-5_dl, 1.73976481843332439e-6_dl, -1.66578804580377724e-6_dl, -8.83571732180485515e-7_dl, -1.15092657356077256e-7_dl, 1.54021206255242797e-7_dl /)
    if(lnam .le. -2.49d0)then
       corr = (5.d0/7.d0/const_pi2)*exp(2.d0*lnam)
       lnp = log((1.d0 - corr)/3.d0)
       return
    endif
    if(lnam .ge. 3.99d0)then
       invam2 = exp(-2.d0*lnam)
       lnp = log(((900.d0/7.d0/const_pi4*const_Riemannzeta5) + ((-2025.d0*const_Riemannzeta7/const_pi4) +(172125.d0*const_RiemannZeta9/2.d0/const_pi4)*invam2)*invam2)) - lnam
       return
    endif
    call get_cheb_value(cp, 20, (lnam - 0.75d0)/(3.25d0), lnp)
  end subroutine fermion_get_lnp

  !!compute d ln rho / d ln am
  subroutine fermion_get_dlnrho(lnam, dlnrho)
    real(dl) lnam, dlnrho, corr, invam2
    real(dl),dimension(20),parameter::cdrho = (/ 0.45853294748427892_dl, 0.59814890542859578_dl, 6.46604723315440388e-2_dl, -0.13048676489827551_dl, -3.72780683875114996e-2_dl, 3.81270384863884870e-2_dl, 1.77677707709482345e-2_dl, -1.05345886610559315e-2_dl, -7.42005947098849181e-3_dl, 2.52170510469850589e-3_dl, 2.80538890088544331e-3_dl, -4.47876939099481114e-4_dl, -9.75418076177981598e-4_dl, 1.33083589852769130e-5_dl, 3.17887727370739473e-4_dl, 3.97485986379130440e-5_dl, -9.38232585120578154e-5_dl, -2.40406623904638736e-5_dl, 2.99526741850547889e-5_dl, 1.34980095732365434e-5_dl /)
    if(lnam .le. -2.49_dl)then
       corr = (5.d0/7.d0/const_pi2)*exp(2.d0*lnam)
       dlnrho = 2.d0*corr/(1.d0 + corr)
       return
    endif
    if(lnam .gt. 3.99d0)then
       invam2 = exp(-2.d0*lnam)
       dlnrho = 1.d0 + ( (-15.d0*const_Riemannzeta5/const_Riemannzeta3) &
            + (45.d0*(10.d0*const_Riemannzeta5**2 + 21.d0*const_Riemannzeta3*const_Riemannzeta7)/4.d0/const_Riemannzeta3**2)*invam2)*invam2
       return
    endif
    call get_cheb_value(cdrho, 20, (lnam-0.75d0)/3.25d0, dlnrho)
  end subroutine fermion_get_dlnrho


    !!compute d ln rho / d ln am
  subroutine fermion_get_dlnp(lnam, dlnp)
    real(dl) lnam, dlnp, corr, invam2
    real(dl), dimension(20),parameter::cdrho = (/ -0.42955931767638317d+00, -0.57929706031950723d+00, -0.10061753508327952d+00, 0.10137410858790824d+00, 0.47912127913089325d-1 , -0.18296247083641887d-1 , -0.17009826170227303d-1 , 0.13239071953237954d-2 , 0.46372977504751667d-2 , 0.67998892886111039d-3 , -0.92969391416750698d-3 , -0.36634500130460678d-3 , 0.10647905926389047d-3 , 0.10129351444516322d-3 , 0.85424263264848575d-5 , -0.16063650866938333d-4 , -0.75048662407985448d-5 , 0.98439667538101790d-7 , 0.17632977601903781d-5 , 0.97383174483000743d-6  /)
    if(lnam .le. -2.49_dl)then
       corr = (5.d0/7.d0/const_pi**2)*exp(2.d0*lnam)
       dlnp = -2.* corr/(1.-corr)
       return
    endif
    if(lnam .gt. 3.99d0)then
       invam2 = exp(-2.d0*lnam)
       dlnp = - ((-2.*2025.d0*const_Riemannzeta7/const_pi4) +(4.*172125.d0*const_RiemannZeta9/2.d0/const_pi4)*invam2) *invam2 &
            / (((900.d0/7.d0/const_pi4*const_Riemannzeta5) + ((-2025.d0*const_Riemannzeta7/const_pi4) +(172125.d0*const_RiemannZeta9/2.d0/const_pi4)*invam2)*invam2)) &
            - 1.d0
       return
    endif
    call get_cheb_value(cdrho, 20, (lnam-0.75d0)/3.25_dl, dlnp)
  end subroutine fermion_get_dlnp
!!------------------------ Boson -----------------------------------------

!!boson energy density and pressure, the unit is the energy density of a massless fermion with the same termperature rho_massless = pi^4/15 T^4
!!input lnam is ln(mc^2/(kT))
!!output lnrho is ln (rho/rho_massless), ln p is ln(p/rho_massless)
!!in the limit lnam -> - infinity, apparently you get lnrho = 0 and lnp = ln(1/3)
  subroutine boson_get_lnrho(lnam, lnrho)
    real(dl) invam2, corr, lnam, lnrho
    real(dl),dimension(20),parameter::crho = (/ 0.85920499059112343_dl, 1.3962721264359432_dl, 0.70432294674814611_dl, 0.14297095706190513_dl, -6.55811511010368081e-2_dl, -4.63591428381798787e-2_dl, 5.01875599077331383e-3_dl, 1.41690006104390458e-2_dl, 2.21177480941542586e-3_dl, -3.70703629807588638e-3_dl, -1.66893509068491972e-3_dl, 7.15614414586352820e-4_dl, 7.23143566964830971e-4_dl, -3.35664246204683407e-5_dl, -2.42009019038368436e-4_dl, -5.66482185018125605e-5_dl, 6.35841417332966958e-5_dl, 3.44090360869310461e-5_dl, -1.09428321341706716e-5_dl, -1.68472313639487885e-5_dl /)
    if(lnam .le. -3.99d0)then
       corr = (5.d0/4.d0/const_pi2)*exp(2.d0*lnam)
       lnrho = log(1.d0 + corr)
       return
    endif
    if(lnam .ge. 3.99d0)then
       invam2 = exp(-2.d0*lnam)
       lnrho = lnam + log((const_Riemannzeta3 * 30.d0/const_pi4) +( (180.d0*const_Riemannzeta5/const_pi4) + (-1350.d0/const_pi4 * const_Riemannzeta7 + (37800.d0 * const_Riemannzeta9/ const_pi4)*invam2)*invam2)*invam2)
       return
    endif
    call get_cheb_value(crho, 20,  lnam/4.d0, lnrho)
  end subroutine boson_get_lnrho

  subroutine boson_get_lnp(lnam, lnp)
    real(dl) invam2, corr, lnam, lnp
    real(dl),dimension(20),parameter::cp = (/ -1.8285425606893564_dl, -1.1995994663208203_dl, -0.63579323906066232_dl, -0.16456435145041079_dl, 3.22461274266371070e-2_dl, 4.02357879959531692e-2_dl, 5.56128174481804183e-3_dl, -7.42808998219154047e-3_dl, -3.68630245376660396e-3_dl, 5.37511637030240472e-4_dl, 1.03818428724796819e-3_dl, 2.31230717643599081e-4_dl, -1.53291700135450192e-4_dl, -1.07074751727548019e-4_dl, -6.96618941049786799e-6_dl, 2.02895088393592875e-5_dl, 1.00029976391815033e-5_dl, -7.29511466983895998e-8_dl, -2.23609786187199272e-6_dl, -1.23155707342790243e-6_dl /)
    if(lnam .le. -3.99d0)then
       corr = (5.d0/4.d0/const_pi2)*exp(2.d0*lnam)
       lnp = log((1.d0 - corr)/3.d0)
       return
    endif
    if(lnam .ge. 3.99d0)then
       invam2 = exp(-2.d0*lnam)
       lnp = log(((120.d0/const_pi4*const_Riemannzeta5) + ((-1800.d0*const_Riemannzeta7/const_pi4) +(75600.d0*const_RiemannZeta9/const_pi4)*invam2)*invam2)) - lnam       
       return
    endif
    call get_cheb_value(cp, 20, lnam/4.d0, lnp)
  end subroutine boson_get_lnp
  

  !!compute d ln rho / d ln am
  subroutine boson_get_dlnrho(lnam, dlnrho)
    real(dl) lnam, dlnrho, corr, invam2
    real(dl),dimension(20),parameter::cdrho = (/ 0.41655827483174224_dl, 0.59177723567654272_dl, 0.13498019867091371_dl, -0.11254573750954061_dl, -7.94769037761112740e-2_dl, 1.86165504389742276e-2_dl, 3.64196607223354896e-2_dl, 3.56021952692404855e-3_dl, -1.31738437483698515e-2_dl, -5.28694710409113346e-3_dl, 3.50481795939059987e-3_dl, 3.05761662838206988e-3_dl, -4.35482712001917550e-4_dl, -1.28139839893068341e-3_dl, -2.23845853250542184e-4_dl, 4.12439057719122951e-4_dl, 1.91001532358664038e-4_dl, -9.65718422490377602e-5_dl, -1.18271863089807312e-4_dl, 1.35735494059749740e-6_dl /)
    if(lnam .le. -3.99_dl)then
       corr = (5.d0/4.d0/const_pi2)*exp(2.d0*lnam)
       dlnrho = 2.d0*corr/(1.d0 + corr)
       return
    endif
    if(lnam .gt. 3.99d0)then
       invam2 = exp(-2.d0*lnam)
       dlnrho = 1.d0 + ( (-12.d0*const_Riemannzeta5/const_Riemannzeta3) &
            + ((72.d0*const_Riemannzeta5**2 + 180.d0*const_Riemannzeta3*const_Riemannzeta7)/const_Riemannzeta3**2)*invam2)*invam2
       return
    endif
    call get_cheb_value(cdrho, 20, lnam/4._dl, dlnrho)
  end subroutine boson_get_dlnrho

  subroutine boson_get_dlnp(lnam, dlnp)
    real(dl) lnam, dlnp, corr, invam2
    real(dl), dimension(20),parameter::cdrho =   (/ -0.38445753351392464d+00, -0.56507922513479980d+00, -0.16911534555220517d+00, 0.70714005528823232d-1 , 0.77731140160035683d-1 , 0.62217410104989388d-2 , -0.22858414171934623d-1 , -0.10462112478909665d-1 , 0.31397644217964569d-2 , 0.42830853648848930d-2 , 0.72075785662759496d-3 , -0.90785063326038519d-3 , -0.55131656392327173d-3 , 0.11882642415894909d-4 , 0.14421947573201171d-3 , 0.60608967122437350d-4 , -0.86129681024491393d-5 , -0.19507626000606408d-4 , -0.90584422331560452d-5 , 0.43136694416906456d-6  /)
    if(lnam .le. -3.99_dl)then
       corr = (5.d0/4.d0/const_pi2)*exp(2.d0*lnam)
       dlnp = -2.* corr/(1.-corr)
       return
    endif
    if(lnam .gt. 3.99d0)then
       invam2 = exp(-2.d0*lnam)
       dlnp = -  (((-2.*1800.d0*const_Riemannzeta7/const_pi4) +(4.*75600.d0*const_RiemannZeta9/const_pi4)*invam2)*invam2)  &
            / ((120.d0/const_pi4*const_Riemannzeta5) + ((-1800.d0*const_Riemannzeta7/const_pi4) +(75600.d0*const_RiemannZeta9/const_pi4)*invam2)*invam2) &
            - 1.d0
       return
    endif
    call get_cheb_value(cdrho, 20, lnam/4._dl, dlnp)
  end subroutine boson_get_dlnp


end module particle_utils


