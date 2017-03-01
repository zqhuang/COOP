program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_REAL,parameter::Tstd = 273.16
  COOP_REAL,parameter::p1 = 500.d0,p2 = 200.d0, p3 = 100.d0
  COOP_REAL,parameter::dP1 = 234.d0, dp2 = 93.4d0, dp3 = 46.68d0
  COOP_REAL::dT_fid_ini = 127.4237
  COOP_INT::i
  COOP_REAL:: b, dT, dT_fid, dTave
  read(*,*) b
  dT_fid = dT_fid_ini
  do i = 1, 5
     dTave = (dT12(b)+ dT23(b)+ dT13(b))/3.d0
     write(*,*) i, sqrt(((dT12(b)-dTave)**2 + (dT23(b)-dTave)**2 + (dT13(b)-dTave)**2)/3.d0), dTave+273.16, (Tstd*dp1/dTave - p1)*(b+dT_fid/dp1)**2,  (Tstd*dp2/dTave - p2)*(b+dT_fid/dp2)**2,  (Tstd*dp3/dTave - p3)*(b+dT_fid/dp3)**2
     dT_fid = dTave 
  enddo

contains

  function dT12( b)
    COOP_REAL:: b, dT12
    dT12  = Tstd*((b+ dT_fid/dp1)**2*dp1 - (b+dT_fid/dp2)**2*dp2)/((b+ dT_fid/dp1)**2*p1 - (b+dT_fid/dp2)**2*p2)
  end function dT12

  function dT23( b)
    COOP_REAL:: b, dT23
    dT23  = Tstd*((b+ dT_fid/dp3)**2*dp3 - (b+dT_fid/dp2)**2*dp2)/((b+ dT_fid/dp3)**2*p3 - (b+dT_fid/dp2)**2*p2)
  end function dT23

  function dT13( b)
    COOP_REAL:: b, dT13
    dT13  = Tstd*((b+ dT_fid/dp3)**2*dp3 - (b+dT_fid/dp1)**2*dp1)/((b+ dT_fid/dp3)**2*p3 - (b+dT_fid/dp1)**2*p1)
  end function dT13


end program Test  
