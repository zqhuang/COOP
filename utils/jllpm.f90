module jllpm_utils
  use sphericalbessel_utils
  implicit none
#include "utils.h"


contains

  function Jllpm(l,lp,m,x)  !!always assume l>=2;x>=0
    !!so far support m<=2 and lp<=6
    integer l,lp,m
    real(dl) lsq, xsq
    real(dl) x, Jllpm
    if(lp.eq.0 .and. m.eq.0)then
       jllpm = sphericalbesselJ(l, x)
       return
    endif
    if(lp .gt. 8 .or. m.gt.2 .or. (lp.gt.6 .and. m.gt.0))then
       jllpm = 0
       return
    endif
    xsq = x*x
    lsq = l*l
    if(x .lt. 1.d-5)then
       select case(m)
       case(0)
          select case(l-lp)
          case(0)
             jllpm = (1.d0 + xsq/(-3.d0 + 4.d0* (lsq+l))*((0.5d0 - lsq-l) &
                  + xsq*( &
                  (((3.d0/32.d0) + (l+lsq)*(-4.d0 + l + lsq)/16.d0))/((l-1.5d0)*(l+2.5d0))) &
                  ))/(1.d0 + 2.d0* l)
          case(1)
             jllpm = x*l/(4.d0*lsq-1.d0)*(1.d0 + xsq/(8.d0* lsq-18.d0)*(2.d0 - lsq &
                  + xsq *(23.d0+2.d0*lsq*(lsq-8.d0))/(48.d0*lsq-300.d0) ))
          case(-1)
             jllpm = -x*(l+1.d0)/(4.d0*lsq+8.d0*l+3.d0)*(1.d0 + xsq/(8.d0*lsq+16.d0*l-10.d0)*(1.d0-lsq -2.d0*l + xsq*(((9.d0 + 2.d0* l*(2.d0 + l)* (-6.d0 + l*(2.d0 + l))))/(12.d0* (-3.d0 + 2.d0* l)*(7.d0 + 2.d0*l)))))
          case(2)
             jllpm = xsq*(l-1.d0)*l/((2.d0 + 4.d0*l)*(3.d0 + 4.d0*(l-2.d0 )*l))*( &
                  1.d0-xsq*(7.d0+2.d0*(l-lsq))/(90.d0+24.d0*(l-lsq)))
          case(-2)
             jllpm = xsq*(2.d0 + (3.d0+ l)*l)/(30.d0 + l*(92.d0 + l*(72.d0 + l*16.d0)))*(1.d0+xsq*((1.d0/8.d0)-l*(0.25d0+l/12.d0))/(l-0.5d0)/(l+3.5d0))
          case(3)
             jllpm = x*xsq*(l-2.d0)*(l-1.d0)*l/(-90.d0+l*(96.d0+l*(336.d0+l*(-384.d0+l*96.d0))))*( &
                  1.d0-xsq*(5.d0+2.d0*l-lsq)/(84.d0+32.d0*l-16.d0*lsq))
          case(-3)
             jllpm = -x*xsq*(1.d0 + l)*(2.d0 + l)*(3.d0 + l)/(630.d0+l*(2112.d0+l*(2064.d0+l*(768.d0+l*96.d0))))*(1.d0+xsq*(2.d0-4.d0*l-lsq)/16.d0/(l-0.5d0)/(l+4.5d0)*(1.d0-xsq*(15.d0 + l*(4.d0 + l)*(-10.d0 + l*(4.d0 + l)))/(40.d0*(l-1.5d0)*(l+5.5d0)* (-2.d0 + l*(4.d0 + l)))))
          case(4)
             jllpm = xsq*xsq* l*(l-1.d0)*(l-2)*(l-3)/(24.d0+48.d0*l)/(105.d0+l*(-352.d0 + l*(344.d0+l*(-128.d0 + l*16.d0)))) * (1.d0 +xsq*(6.5d0+3.d0*l-lsq)/(20.d0*(l-4.5d0)*(l+1.5d0)))
          case(-4)
             jllpm = xsq*xsq*(l+1.d0)*(l+2.d0)*(l+3.d0)*(l+4.d0)/(24.d0*32.d0)/((l+0.5d0)*(l+1.5d0)*(l+2.5d0)*(l+3.5d0)*(l+4.5d0))*(1.d0+  xsq*(2.5d0 - 5.d0* l - lsq)/(20.d0* (l-0.5d0)*(l+5.5d0)))
          case(5)
             jllpm = xsq*xsq*x*l*(l-1.d0)*(l-2.d0)*(l-3.d0)*(l-4.d0)/(120.d0*64.d0)/((l-4.5d0)*(l-3.5d0)*(l-2.5d0)*(l-1.5d0)*(l-0.5d0)*(l+0.5d0))
          case(-5)
             jllpm = -xsq*xsq*x*(l+1.d0)*(l+2.d0)*(l+3.d0)*(l+4.d0)*(l+5.d0)/(120.d0*64.d0)/((l+0.5d0)*(l+1.5d0)*(l+2.5d0)*(l+3.5d0)*(l+4.5d0)*(l+5.5d0))
          case(6)
             jllpm = xsq*xsq*xsq*l*(l-1.d0)*(l-2.d0)*(l-3.d0)*(l-4.d0)*(l-5.d0)/(720.d0*128.d0)/((l-5.5d0)*(l-4.5d0)*(l-3.5d0)*(l-2.5d0)*(l-1.5d0)*(l-0.5d0)*(l+0.5d0))
          case(-6)
             jllpm = xsq*xsq*xsq*(l+1.d0)*(l+2.d0)*(l+3.d0)*(l+4.d0)*(l+5.d0)*(l+6.d0)/(720.d0*128.d0)/((l+0.5d0)*(l+1.5d0)*(l+2.5d0)*(l+3.5d0)*(l+4.5d0)*(l+5.5d0)*(l+6.5d0))
          case default
             jllpm = 0.d0
          end select
       case(1)
          select case(l-lp)
          case(0)
             jllpm = (1.d0+xsq/(4.d0*(l+lsq)-3.d0)*(1.5d0-l-lsq + xsq*(15.d0+2.d0*(l-2.d0)*(lsq+l)*(l+3.d0)/32.d0/(l-1.5d0)/(l+2.5d0))))/(1.d0+2.d0*l)
          case(1)
             jllpm = x*sqrt_table(l-1)*sqrt_table(l+1)/(4.d0*lsq-1.d0)*(1.d0+xsq/(8.d0*lsq-18.d0)*(3.d0-lsq+xsq*(45.d0+lsq*(-20.d0+lsq*2.d0))/(-300.d0+48.d0*lsq)))
          case(-1)
             jllpm = -x*sqrt_table(l)*sqrt_table(l+2)/(4.d0*lsq+8.d0*l+3.d0)*(1.d0+xsq/(4.d0*lsq+8.d0*l-5.d0)*(-lsq/2.d0-l+1.d0+xsq*(27.d0+2.d0*l*(lsq-4.d0)*(l+4.d0))/96.d0/(l-1.5d0)/(l+3.5d0)))
          case(2)
             jllpm = xsq* sqrt_table(l-2)*sqrt_table(l-1)*sqrt_table(l)*sqrt_table(l+1)/(6.d0 + l*(- 4.d0 + l*(- 24.d0+l*16.d0)))*(1.d0-xsq*(9.d0+2.d0*(l-lsq))/(90.d0+24.d0*(l-lsq)))
          case(-2)
             jllpm = xsq* sqrt_table(l+3)*sqrt_table(l)*sqrt_table(l+1)*sqrt_table(l+2)/(30.d0+l*(92.d0+l*(72.d0+l*16.d0)))*(1.d0+xsq*(5.d0-2.d0*l*(l+3.d0))/(24.d0*(l-0.5d0)*(l+3.5d0)))
          case(3)
             jllpm = x*xsq*(l-1.d0)*sqrt_table(l)*sqrt_table(l+1)*sqrt_table(l-2)*sqrt_table(l-3)/(-90.d0+l*(96.d0+l*(336.d0+l*(-384.d0+l*96.d0))))*(1.d0-xsq*(6.d0+2.d0*l-lsq)/(84.d0+32.d0*l-16.d0*lsq))
          case(-3)
             jllpm =-x*xsq*(l+2.d0)*sqrt_table(l)*sqrt_table(l+1)*sqrt_table(l+3)*sqrt_table(l+4)/(630.d0+l*(2112.d0+l*(2064.d0+l*(768.d0+l*96.d0))))*(1.d0+xsq*(3.d0-l*(l+4.d0))/(16.d0*(l-0.5d0)*(l+4.5d0)))
          case(4)
             jllpm = ((-2.d0 + l)*(-1.d0 + l)*sqrt_table(l)*sqrt_table(l + 1)*sqrt_table(l-4)*sqrt_table(l-3)*xsq*xsq)/(24.d0*32.d0)/((-3.5d0 + l)*(-2.5d0 + l)*(-1.5d0 +l)*(-0.5d0 + l)*(0.5d0 + l))
          case(-4)
             jllpm = ((2.d0 + l)*(3.d0 + l)*sqrt_table(l)*sqrt_table(l+1)*sqrt_table(l+4)*sqrt_table(l+5)*xsq*xsq)/(24.d0*32.d0)/((0.5d0 + l)*(1.5d0 + l)*(2.5d0 + l)*(3.5d0 + l)*(4.5d0 + l))
          case default
             jllpm = 0.d0
          end select
       case(2)
          select case(l-lp)
          case(0)
             jllpm = (1.d0+xsq/(4.d0*(l+lsq)-3.d0)*(4.5d0-l-lsq+xsq*(75.d0+2.d0*(l-3.d0)*l*(l+1.d0)*(l+4.d0))/32.d0/(l-1.5d0)/(l+2.5d0)))/(1.d0+2.d0*l)
          case(1)
             jllpm = x*sqrt_table(l-2)*sqrt_table(l+2)/(4.d0*lsq-1.d0)*(1.d0+xsq/(8.d0*lsq-18.d0)*(6.d0-lsq+xsq*(135.d0+lsq*(-32.d0+lsq*2.d0))/(48.d0*lsq-300.d0)))
          case(-1)
             jllpm = -x*sqrt_table(l+3)*sqrt_table(l-1)/(4.d0*lsq+8.d0*l+3.d0)*(1.d0+xsq/(4.d0*lsq+8.d0*l-5.d0)*(-lsq/2.d0-l+2.5d0+xsq*(105.d0+2.d0*l*(l+2.d0)*(-14.d0+l*(2.d0+l)))/96.d0/(l-1.5d0)/(l+3.5d0)))
          case(2)
             jllpm = xsq*sqrt_table(l+2)*sqrt_table(l+1)*sqrt_table(l-2)*sqrt_table(l-3)/(6.d0 + l*(- 4.d0 + l*(- 24.d0+l*16.d0)))*(1.d0-xsq*(15.d0+2.d0*(l-lsq))/(90.d0+24.d0*(l-lsq)))
          case(-2)
             jllpm = xsq* sqrt_table(l+3)*sqrt_table(l+4)*sqrt_table(l)*sqrt_table(l-1)/(30.d0+l*(92.d0+l*(72.d0+l*16.d0)))*(1.d0+xsq*(11.d0-2.d0*l*(l+3.d0))/24.d0/(l-0.5d0)/(l+3.5d0))
          case(3)
             jllpm = x*xsq*sqrt_table(l+2)*sqrt_table(l+1)*sqrt_table(l)*sqrt_table(l-2)*sqrt_table(l-3)*sqrt_table(l-4)/(-90.d0+l*(96.d0+l*(336.d0+l*(-384.d0+l*96.d0))))*(1.d0+((9.d0 - (-2.d0 + l)*l)*xsq)/(16.d0*(-3.5d0 + l)*(1.5d0 + l)))
          case(-3)
             jllpm = -x*xsq*sqrt_table(l+5)*sqrt_table(l+4)*sqrt_table(l+3)*sqrt_table(l+1)*sqrt_table(l-1)*sqrt_table(l)/(630.d0+l*(2112.d0+l*(2064.d0+l*(768.d0+l*96.d0))))*(1.d0-((-6.d0 + l*(4.d0 + l))*xsq)/(16.*(-0.5d0 + l)*(4.5d0 + l)))
          case(4)
             jllpm = xsq*xsq*sqrt_table(l-5)*sqrt_table(l-4)*sqrt_table(l-3)*sqrt_table(l-2)*sqrt_table(l-1)*sqrt_table(l)*sqrt_table(l+1)*sqrt_table(l+2)/(24.d0*32.d0)/((l-3.5d0)*(l-2.5d0)*(l-1.5d0)*(l-0.5d0)*(l+0.5d0))
          case(-4)
             jllpm = xsq*xsq*sqrt_table(l-1)*sqrt_table(l+1)*sqrt_table(l)*sqrt_table(l+2)*sqrt_table(l+3)*sqrt_table(l+4)*sqrt_table(l+5)*sqrt_table(l+6)/(24.d0*32.d0)/((l+0.5d0)*(l+1.5d0)*(l+2.5d0)*(l+3.5d0)*(l+4.5d0))
          case default
             jllpm = 0.d0
          end select
       end select
       return
    endif
    if(xsq .lt. l)then
       select case(m)
       case(0)
          select case(lp)
          case(0)
             jllpm =  dexp(-lngamma_hi_table(l+2)+l*(dlog(x)-const_ln2)+dlog(const_sqrtpi/2.d0))  * ( &
                  1.d0  &
                  - xsq/(l+1.5d0)/4.d0*(1.d0 &
                  - xsq/(l+2.5d0)/8.d0*(1.d0 &
                  - xsq/(l+3.5d0)/12.d0*(1.d0 &
                  - xsq/(l+4.5d0)/16.d0*(1.d0 &
                  - xsq/(l+5.5d0)/20.d0 * (1.d0 &
                  ))))))
          case(1)
             jllpm = dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi/4.d0)) * ( &
                  l &
                  - xsq / (l + 1.5d0) / 4.d0 *  ( l + 2 &
                  - xsq / (l + 2.5d0) / 8.d0 *  ( l + 4 &
                  - xsq / (l + 3.5d0) / 12.d0 * ( l + 6 &
                  - xsq / (l + 4.5d0) / 16.d0 * ( l + 8 &
                  - xsq / (l + 5.5d0) / 20.d0 * ( l + 10 &
                  ))))))
          case(2)
             jllpm =  dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi*3.d0/16.d0)) * ( &
                  l*(l-1.d0) &
                  - xsq / ( l + 1.5d0) / 4.d0 * ( l*(l+5.d0/3.d0) &
                  - xsq / ( l + 2.5d0) / 8.d0 * ( 16.d0/3.d0 + l*(13.d0/3.d0+l) &
                  - xsq / ( l + 3.5d0) / 12.d0 * (16.d0 + l * (l+7.d0) &
                  - xsq / ( l + 4.5d0) / 16.d0 * (32.d0 + l * (l+29.d0/3.d0) &
                  - xsq / ( l + 5.5d0) / 20.d0 * (160.d0/3.d0 + l * (l+37.d0/3.d0) &
                  ))))))
          case(3)
             jllpm =  dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi*5.d0/32.d0)) * ( &
                  l*(l-1.d0)*(l-2) &
                  - xsq/(l+1.5d0)/4.d0 * ( l*(l-1)*(l+1.6d0) &
                  - xsq/(l+2.5d0)/8.d0 * ( l*(l+2)*(l+2.2d0) &
                  - xsq/(l+3.5d0)/12.d0*( (l+4)*(4.8d0+l*(l+3.8d0)) &
                  - xsq/(l+4.5d0)/16.d0*((l+6)*(12.8d0+l*(l+5.4d0)) &
                  - xsq/(l+5.5d0)/20.d0*((l+8)*(24.d0+l*(l+7.d0)) &
                  ))))))
          case(4)
             jllpm =  dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi*35.d0/256.d0)) * ( &
                  l*(l-1.d0)*(l-2)*(l-3) &
                  - xsq/(l+1.5d0)/4.d0 * ( l*(l-1.d0)*(l-2)*(l+11.d0/7.d0) &
                  - xsq/(l+2.5d0)/8.d0 * ( l*(l-1)* ( 146.d0/35.d0 + (l+29.d0/7.d0)*l) &
                  - xsq/(l+3.5d0)/12.d0 * (l*( 558.d0/35.d0 + l*(673.d0/35.d0 + l*(54.d0/7.d0+l))) &
                  - xsq/(l+4.5d0)/16.d0*( 3072.d0/35.d0 + l*(3758.d0/35.d0+l*(1921.d0/35.d0+l*(86.d0/7.d0+l))) &
                  - xsq/(l+5.5d0)/20.d0*(3072.d0/7.d0+l*(314.d0+l*(107.d0+l*(118.d0/7.d0+l))) &
                  ))))))
          case(5)
             jllpm =  dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi*63.d0/512.d0)) * ( &
                  l*(l-1.d0)*(l-4)*(l-3)*(l-2) &
                  - xsq/(l+1.5d0)/4.d0 * ( (l+14.d0/9.d0)*(l-3)*(l-2)*(l-1)*l &
                  - xsq/(l+2.5d0)/8.d0 * ( (256.d0/63.d0+l*(37.d0/9.d0+l))*(l-2)*(l-1)*l &
                  - xsq/(l+3.5d0)/12.d0 * ( (158.d0/21.d0+l*(17.d0/3.d0+l))*(l+2)*(l-1)*l &
                  - xsq/(l+4.5d0)/16.d0 * ( l*(l+4)*(1158.d0/63.d0 + l*(1333.d0/63.d0 + l*(74.d0/9.d0+l))) &
                  - xsq/(l+5.5d0)/20.d0 * ( (l+6)*(5120.d0/63.d0+l*(6182.d0/63.d0+l*(3253.d0/63.d0+l*(106.d0/9.d0+l))))&
                  - xsq/(l+6.5d0)/24.d0 * ( (l+8)*(7680.d0/21.d0+l*(5378.d0/21.d0+l*(1959.d0/21.d0+l*(46.d0/3.d0+l))))&
                  )))))))
          case(6)
             jllpm =   dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi*231.d0/2048.d0)) * ( &
                  l*(l-1.d0)*(l-2)*(l-3)*(l-4)*(l-5) &
                  - xsq/(l+1.5d0)/4.d0 * ( l*(l-1.d0)*(l-4)*(l-3)*(l-2)*(l+17.d0/11.d0) &
                  - xsq/(l+2.5d0)/8.d0 * ( l*(l-1.d0)*(l-2)*(l-3)*(4.d0+l*(45.d0/11.d0+l)) & 
                  - xsq/(l+3.5d0)/12.d0 * ( l*(l-1.d0)*(l-2)*(1124.d0/77.d0+l*(205.d0/11.d0+l*(84.d0/11.d0+l))) &
                  - xsq/(l+4.5d0)/16.d0 * ( l*(l-1)*(5336.d0/77.d0+l*(7758.d0/77.d0+l*(589.d0/11.d0+l*(134.d0/11.d0+l)))) & 
                  - xsq/(l+5.5d0)/20.d0 * ( l*(31720.d0/77.d0+l*(49498.d0/77.d0+l*(30875.d0/77.d0+l*(1335.d0/11.d0+l*(195.d0/11.d0+l))))) &
                  - xsq/(l+6.5d0)/24.d0 * (245760.d0/77.d0 + l*(369416.d0/77.d0+l*(255098.d0/77.d0+l*(92515.d0/77.d0+l*(2615.d0/11.d0+l*(267.d0/11.d0+l))))) &
                  )))))))
          case(7)
             jllpm = dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi*429.d0/4096.d0))*( &
                  l*(l-1.d0)*(l-2)*(l-3)*(l-4)*(l-5)*(l-6) &
                  -xsq/(l+1.5d0)/4.d0*(l*(l-1.d0)*(-5+l)*(-4+l)*(-3+l)*(-2+l)*(l + 20.d0/13.d0) &
                  -xsq/(l+2.5d0)/8.d0*( l*(l-1.d0)*(l-2)*(l-3)*(l-4)*(l*(l+53.d0/13.d0)+566.d0/143.d0) &
                  - xsq/(l+3.5d0)/12.d0*( l*(l-1.d0)*(l-2)*(l-3)*(l*(l*(l+99.d0/13.d0) + 2644.d0/143.d0 )+2048.d0/143.d0) &
                  -xsq/(l+4.5d0)/16.d0*( l*(l-1.d0)*(l-2)*(l+2)*(l*(l*(l+132.d0/13.d0)+4705.d0/143.d0)+4796.d0/143.d0) &
                  -xsq/(l+5.5d0)/20.d0*( l*(l-1)*(l+4.d0)*(l*(l*(l*(l+178.d0/13.d0)+9433.d0/143.d0)+18978.d0/143.d0)+13880.d0/143.d0) &
                  - xsq/(l+6.5d0)/24.d0*(l*(l*(l+6.d0)*(l*(l*(l*(l+237.d0/13.d0)+18203.d0/143.d0)+61057.d0/143.d0)+98990.d0/143.d0)+64760.d0/143.d0) &
                  )))))))
          case(8)
             jllpm = dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi*6435.d0/65536.d0))*( &
                  l*(l-1.d0)*(l-2)*(l-3)*(l-4)*(l-5)*(l-6)*(l-7) &
                  - xsq/(l+1.5d0)/4.d0*( l*(l-1.d0)*(l-2)*(l-3)*(l-4)*(l-5)*(l-6)*(l+23.d0/15.d0) &
                  - xsq/(l+2.5d0)/8.d0*( l*(l-1.d0)*(l-2)*(l-3)*(l-4)*(l-5)*(l*(l+61.d0/15.d0)+766.d0/195.d0) &
                  - xsq/(l+3.5d0)/12.d0*( l*(l-1.d0)*(l-2)*(l-3)*(l-4)*(10106.d0/715.d0 + l*(1195.d0/65.d0+l*(38.d0/5.d0+l))) &
                  - xsq/(l+4.5d0)/16.d0*( l*(l-1.d0)*(l-2)*(l-3)*(140872.d0/2145.d0 + l*(211010.d0/2145.d0+l*(10329.d0/195.d0+l*(182.d0/15.d0+l)))) &
                  - xsq/(l+5.5d0)/20.d0*( l*(l-1.d0)*(l-2)*(161032/429.d0+l*(265242/429.d0+l*(168821/429.d0+l*(4691.d0/39.d0+l*(53.d0/3.d0+l)) ))) &
                  - xsq/(l+6.5d0)/24.d0*( l*(l-1.d0)*(1832080.d0/715.d0+l*(3215372.d0/715.d0+l*(2298960.d0/715.d0+l*(846005.d0/715.d0+l*(15335.d0/65.d0+l*(121.d0/5.d0+l)))))) &
                  - xsq/(l+7.5d0)/28.d0*( l* (44239440./2145.d0 + l* (80163724./2145.d0 + l*(62311732./2145.d0 + l * (26154065./2145.d0 + l*(6364120./2145.d0 + l*(81606/195.d0 + l * (476/15.d0 + l)))))))  &
                  ))))))))
          case default
             jllpm = 0
          end select
       case(1)
          select case(lp)
          case(1)
             jllpm =  dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi/4.d0/const_sqrt2)) * sqrt_table(l)*sqrt_table(l+1) * ( &
                  1.d0  &
                  - xsq/(l+1.5d0)/4.d0*(1.d0 &
                  - xsq/(l+2.5d0)/8.d0*(1.d0 &
                  - xsq/(l+3.5d0)/12.d0*(1.d0 &
                  - xsq/(l+4.5d0)/16.d0*(1.d0 &
                  - xsq/(l+5.5d0)/20.d0 * (1.d0 &
                  ))))))
          case(2)
             jllpm =  dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi/8.d0*const_sqrt3/const_sqrt2)) * sqrt_table(l)*sqrt_table(l+1) * ( &                
                  l-1 &
                  - xsq / (l + 1.5d0) / 4.d0 *  ( l + 1 &
                  - xsq / (l + 2.5d0) / 8.d0 *  ( l + 3 &
                  - xsq / (l + 3.5d0) / 12.d0 * ( l + 5 &
                  - xsq / (l + 4.5d0) / 16.d0 * ( l + 7 &
                  - xsq / (l + 5.5d0) / 20.d0 * ( l + 9 &
                  ))))))
          case(3)
             jllpm =  dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi*5.d0/64.d0*const_sqrt3)) * sqrt_table(l)*sqrt_table(l+1) * ( &     
                  (l-2)*(l-1) &
                  - xsq/(l+1.5d0)/4.d0 * ( (l-1)*(l+1.2d0) &
                  - xsq/(l+2.5d0)/8.d0 * ( 2.d0+l*(3.4d0+l) &
                  - xsq/(l+3.5d0)/12.d0 * (11.6d0+l*(l+6.6d0) &
                  - xsq/(l+4.5d0)/16.d0*( 27.6d0+l*(l+9.8d0) &
                  - xsq/(l+5.5d0)/20.d0 * (50.d0+l*(l+13.d0) &                 
                  ))))))
          case(4)
             jllpm =  dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi*7.d0/128.d0*const_sqrt5)) * sqrt_table(l)*sqrt_table(l+1) * ( &    
                  (l-3)*(l-2)*(l-1.d0) &
                  - xsq/(l+1.5d0)/4.d0*( (l-2)*(l-1)*(l+9.d0/7.d0)  &
                  - xsq/(l+2.5d0)/8.d0*(  (l-1)*(l+1)*(l+18.d0/7.d0) &
                  - xsq/(l+3.5d0)/12.d0*( (l+3)*(2.d0+l*(l+27.d0/7.d0))  &
                  - xsq/(l+4.5d0)/16.d0*( (l+5)*(78.d0/7.d0+l*(43.d0/7.d0+l))  &
                  - xsq/(l+5.5d0)/20.d0*( (l+7)*(174.d0/7.d0+l*(l+59.d0/7.d0))  &
                  ))))))
          case(5)
             jllpm =  dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi*21.d0/512.d0*dsqrt(7.5d0))) * sqrt_table(l)*sqrt_table(l+1) * ( &
                  (l-1.d0)*(l-4)*(l-3)*(l-2) &
                  - xsq/(l+1.5d0)/4.d0*( (l-3)*(l-2)*(l-1.d0)*(l+4.d0/3.d0)  &
                  - xsq/(l+2.5d0)/8.d0*(  (l-2)*(l-1)*(20.d0/7.d0+l*(11.d0/3.d0+l)) &
                  - xsq/(l+3.5d0)/12.d0*( (l-1)*(8.d0+l*(102.d0/7.d0+l*(l+7.d0))) &
                  - xsq/(l+4.5d0)/16.d0*( 168/7.d0+l*(1382.d0/21.d0 + l*(309.d0/7.d0+l*(34.d0/3.d0+l))) &
                  - xsq/(l+5.5d0)/20.d0*( 2248.d0/7.d0+l*(6070.d0/21.d0+l*(725.d0/7.d0+l*(50.d0/3.d0+l))) &
                  ))))))
          case(6)
             jllpm =  dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi*33.d0/1024.d0*dsqrt(10.5d0))) * sqrt_table(l)*sqrt_table(l+1) * ( &
                  (l-1.d0)*(l-2)*(l-3)*(l-4)*(l-5) &
                  - xsq/(l+1.5d0)/4.d0 * ( (l-1.d0)*(l-4)*(l-3)*(l-2)*(l+15.d0/11.d0) &
                  - xsq/(l+2.5d0)/8.d0 * ( (l-1.d0)*(l-2)*(l-3)*(100.d0/33.d0+l*(41.d0/11.d0+l)) &
                  - xsq/(l+3.5d0)/12.d0 * ( (l-1.d0)*(l-2)*(l+1.d0)*(100.d0/11.d0+l*(67.d0/11.d0+l)) &
                  - xsq/(l+4.5d0)/16.d0 * ( (l-1)*(l+3)*(120.d0/11.d0+l*(222.d0/11.d0+l*(93/11.d0+l))) &
                  - xsq/(l+5.5d0)/20.d0 * ( (l+5)*(24.d0 + l*(2350.d0/33.d0 + l*(1555.d0/33.d0+l*(130.d0/11.d0+l)))) & 
                  - xsq/(l+6.5d0)/24.d0 * ( (l+7)*(3400.d0/11.d0+l*(3002.d0/11.d0+l*(1089.d0/11.d0+l*(178.d0/11.d0+l)))) &
                  )))))))
          case default
             jllpm = 0
          end select
       case(2)
          select case(lp)
          case(2)
             jllpm =  dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi*const_sqrt3/const_sqrt2/16.d0)) * sqrt_table(l)*sqrt_table(l+1)*sqrt_table(l+2)*sqrt_table(l-1) * ( &
                  1.d0  &
                  - xsq/(l+1.5d0)/4.d0*(1.d0 &
                  - xsq/(l+2.5d0)/8.d0*(1.d0 &
                  - xsq/(l+3.5d0)/12.d0*(1.d0 &
                  - xsq/(l+4.5d0)/16.d0*(1.d0 &
                  - xsq/(l+5.5d0)/20.d0 * (1.d0 &
                  ))))))
          case(3)
             jllpm =  dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi*dsqrt(7.5d0)/32.d0)) * sqrt_table(l)*sqrt_table(l+1)*sqrt_table(l+2)*sqrt_table(l-1) * ( &
                  l-2 &
                  - xsq/(l+1.5d0)/4.d0*(l &
                  - xsq/(l+2.5d0)/8.d0*(l+2 &
                  - xsq/(l+3.5d0)/12.d0*(l+4 &
                  - xsq/(l+4.5d0)/16.d0*(l+6 &
                  - xsq/(l+5.5d0)/20.d0 * (l+8 &
                  ))))))
          case(4)
             jllpm =  dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi*dsqrt(2.5d0)*7.d0/128.d0)) * sqrt_table(l)*sqrt_table(l+1)*sqrt_table(l+2)*sqrt_table(l-1) * ( &
                  (l-3)*(l-2) &
                  - xsq/(l+1.5d0)/4.d0*((l-2)*(l+3.d0/7.d0) &
                  - xsq/(l+2.5d0)/8.d0*(-6.d0/7.d0+l*(13.d0/7.d0+l) &
                  - xsq/(l+3.5d0)/12.d0*(6.d0 +l*(37.d0/7.d0+l) &
                  - xsq/(l+4.5d0)/16.d0*(138.d0/7.d0+l*(61.d0/7.d0+l) &
                  - xsq/(l+5.5d0)/20.d0 * (282.d0/7.d0+l*(85.d0/7.d0+l) &
                  ))))))
          case(5)
             jllpm =  dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi*dsqrt(52.5d0)*3.d0/256.d0)) * sqrt_table(l)*sqrt_table(l+1)*sqrt_table(l+2)*sqrt_table(l-1) * ( &
                  (l-2.d0)*(l-3)*(l-4) &
                  - xsq/(l+1.5d0)/4.d0*((l-3)*(l-2)*(l+2.d0/3.d0) &
                  - xsq/(l+2.5d0)/8.d0*((l-2)*l*(l+7.d0/3.d0) &
                  - xsq/(l+3.5d0)/12.d0*((l+2)*(-2+l*(l+3)) &
                  - xsq/(l+4.5d0)/16.d0*((l+4)*(6.d0+l*(17.d0/3.d0+l)) &
                  - xsq/(l+5.5d0)/20.d0 * ((l+6)*(58.d0/3.d0+l*(25.d0/3.d0+l))  &
                  ))))))
          case(6)
             jllpm =  dexp(-lngamma_hi_table(l+2)+(l-lp)*(dlog(x)-const_ln2)+dlog(const_sqrtpi*dsqrt(105.d0)*33.d0/4096.d0)) * sqrt_table(l)*sqrt_table(l+1)*sqrt_table(l+2)*sqrt_table(l-1) * ( &
                  (l-2.d0)*(l-3)*(l-4)*(l-5) &
                  - xsq/(l+1.5d0)/4.d0*((l-3.d0)*(l-4)*(l-2)*(l+9.d0/11.d0)  &
                  - xsq/(l+2.5d0)/8.d0*((l-3)*(l-2)*(20.d0/33.d0+l*(29.d0/11.d0+l)) &
                  - xsq/(l+3.5d0)/12.d0*((l-2)*(-20.d0/11.d0+l*(69.d0/11.d0+l*(60.d0/11.d0+l))) &
                  - xsq/(l+4.5d0)/16.d0*(-216.d0/11.d0+l*(98.d0/11.d0+l*(269.d0/11.d0+l*(102.d0/11.d0+l))) &
                  - xsq/(l+5.5d0)/20.d0 * (3960.d0/33.d0+l*(5798.d0/33.d0+l*(2663.d0/33.d0+l*(498.d0/33.d0+l))) &
                  - xsq/(l+6.5d0)/24.d0 * (8744.d0/11.d0+l*(6370.d0/11.d0+l*(1805.d0/11.d0+l*(230.d0/11.d0+l))) &
                  )))))))
          case default
             jllpm = 0
          end select
       end select
       return
    end if


    select case(m)
    case(0)
       select case(lp)
       case(0)
          jllpm = SphericalBesselJ(l,x)
       case(1)
          jllpm = (l* SphericalBesselJ(l-1,x) - (l+1.d0)*SphericalBesselJ(l+1,x))/(2.d0*l+1.d0)
       case(2) 
          jllpm = ( (1.5d0*(lsq-l)/x-x) * SphericalBesselJ(l-1,x) &
               + (1.5d0*(lsq+3.d0*l+2.d0)/x-x) * SphericalBesselJ(l+1,x) &
               )/(2.d0*l+1.d0)
       case(3)
          jllpm = ( (2.5d0*(l-2)*(l-1.d0)*l/xsq - l+5.d0)*SphericalBesselJ(l-1,x) &
               - (2.5d0*(l+1.d0)*(l+2.d0)*(l+3.d0)/xsq-l-6.d0)* SphericalBesselJ(l+1,x) &
               )/(2.d0*l+1.d0)
       case(4)          
          jllpm =  ( (((35.d0/8.d0)*(l-3)*(l-2)*(lsq-l)/xsq - 5.d0*(lsq-l+7.d0))/x + x)*SphericalBesselJ(l-1,x) &
               + (((35.d0/8.d0)*(l+1.d0)*(l+2.d0)*(l+3.d0)*(l+4.d0)/xsq-5.d0*(lsq+3.d0*l+9.d0))/x + x ) *SphericalBesselJ(l+1,x) &
               )/(2.d0*l+1.d0)
       case(5)       
          jllpm = ( (((63.d0/8.d0)*(l-4)*(l-3)*(l-2)*(lsq-l)/xsq-7.d0*(-45.d0+l*(2.d0+l*(l-12.d0))))/xsq + l-14.d0) *SphericalBesselJ(l-1,x) &
               + (((-63.d0/8.d0)*(l+1.d0)*(l+2.d0)*(l+3.d0)*(l+4.d0)*(l+5.d0)/xsq+ 7.d0*(60.d0+l*(29.d0+l*(15.d0+l))))/xsq-l-15.d0) * SphericalBesselJ(l+1,x) &
               )/(2.d0*l+1.d0)
       case(6)
          jllpm = ( ((((231.d0/16.d0)*(l-5)*(l-4)*(l-3)*(l-2)*(lsq-l)/xsq+(-63.d0/8.d0)*(440.d0+l*(26.d0+3.d0*l*(55.d0+l*(-6.d0+l)))))/xsq+10.5d0*(18.d0+lsq-l))/x - x) * SphericalBesselJ(l-1,x) &
               + ((((231.d0/16.d0)*(l+1.d0)*(l+2.d0)*(l+3.d0)*(l+4.d0)*(l+5.d0)*(l+6.d0)/xsq+(-63.d0/8.d0)*(600.d0+l*(370.d0+3.d0*l*(79.d0+l*(10.d0+l)))))/xsq+10.5d0*(20.d0+lsq+l*3.d0))/x - x) * SphericalBesselJ(l+1,x) &
               )/(2.d0*l+1.d0)
       case(7)
          jllpm = ( &
               ((((429.d0/16.d0)*(l-6)*(l-5)*(l-4)*(l-3)*(l-2)*(lsq-l)/xsq+(-99.d0/8.d0)*(-3640.d0 + l*(-578.d0 + l*(-1723.d0 + 3.d0*l*(61.d0 + (-23.d0 + l)*l)))))/xsq+4.5d0*(-616.d0 + 3.d0*l*(2.d0 + (-25.d0 + l)*l)))/xsq+27.d0-l)*SphericalBesselJ(l-1,x) &
               + ((((-429.d0/16.d0)*(l+1.d0)*(l+2.d0)*(l+3.d0)*(l+4.d0)*(l+5.d0)*(l+6.d0)*(l+7.d0)/xsq+(99.d0/8.d0)*(5040.d0 + l*(3708.d0 + l*(2716.d0 + 3.d0* l* (163.d0 + l* (28.d0 + l))))))/xsq-4.5d0*(700.d0 + 3.d0*l*(55.d0 + l* (28.d0 + l))))/xsq+l+28.d0) * SphericalBesselJ(l+1,x) &
               )/(2.d0*l+1.d0)
       case(8) 
          jllpm = ( &
               (((( (6435.d0/128.d0)*(l-7)*(l-6)*(l-5)*(l-4)*(l-3)*(l-2)*(lsq-l)/xsq - (13728.d0/128.d0)*(6300.d0 + l*(1590.d0+ l*(3559.d0 + l*(-315.d0 + l*(220.d0 + (-15.d0 + l)*l))))))/xsq +(3168.d0/128.d0)*(1820.d0 + l*(86.d0 + 3.d0*l*(115.d0 + (-6.d0 + l)*l))))/xsq - (2304.d0/128.d0)*(33.d0 + lsq-l))/x + x)*SphericalBesselJ(l-1,x) &
               +(((( (6435.d0/128.d0)*(1.d0 + l)*(2.d0 + l)*(3.d0 + l)*(4.d0 + l)*(5.d0 + l)*(6.d0 + l)*(7.d0 + l)*(8.d0 + l)/xsq - (13728.d0/128.d0)*(8820.d0 + l*(7434.d0 + l*(5989.d0 + l*(1365.d0 + l*(310.d0 + l*(21.d0 + l)))))))/xsq +(3168.d0/128.d0)*(2100.d0 + l*(670.d0 + 3.d0*l*(139.d0 + l*(10.d0 + l)))))/xsq - (2304.d0/128.d0)*(35.d0 + l*(3.d0 + l)))/x + x)*SphericalBesselJ(l+1,x) &
               )/(2.d0*l+1.d0)
       case(9)
          jllpm = ( &
               ( ((((12155.d0/128.d0)*(-8 + l)*(-7 + l)*(-6 + l)*(-5 + l)*(-4 + l)*(-3 + l)*(-2 + l)*(lsq-l)/xsq - (22880.d0/128.d0)*(-64260.d0 + l*(-21822.d0 + l*(-41867.d0 + l*(2644.d0 + l*(-3455.d0 + l*(277.d0 + (-38.d0 + l)*l)))))))/xsq + (13728.d0/128.d0) * (-7560.d0 + l*(-816.d0 + l*(-1940.d0 + l*(95.d0 + (-40.d0 + l)*l)))))/xsq - (2816.d0/128.d0)*(-585.d0 + l*(2.d0 + (-42.d0 + l)*l)))/xsq + (-44.d0 + l)) * SphericalBesselJ(l-1,x) &
               - ( ((((12155.d0/128.d0)*(1.d0 + l)*(2.d0 + l)*(3.d0 + l)*(4.d0 + l)*(5.d0 + l)*(6.d0 + l)*(7.d0 + l)*(8.d0 + l)*(9.d0 + l)/xsq - (22880.d0/128.d0)*(90720.d0 + l*(85284.d0 + l*(73890.d0 + l*(20029.d0 + l*(5445.d0 + l*(526.d0 + l*(45.d0 + l))))))))/xsq + (13728.d0/128.d0)*(8820.d0 + l*(3514.d0 + l*(2475.d0 + l*(265.d0 + l*(45.d0 + l))))))/xsq - (2816.d0/128.d0)*(630.d0 + l*(89.d0 + l*(45.d0 + l))))/xsq + (45.d0 + l)) * SphericalBesselJ(l+1,x) &
               )/(2.d0*l+1.d0)
       case(10)
          jllpm = ( &
               ((((((46189.d0/256.d0) *(-9 + l)* (-8 + l)* (-7 + l)* (-6 + l)* (-5 + l)* (-4 + l)*(-3 + l)* (-2 + l)* (lsq-l)/xsq - &
               (24310.d0/256.d0)*(2298240.d0 +  l*(967056.d0 + l* (1687484.d0 +  l*(-57452.d0 + 5.d0* l* (34433.d0 + l*(-2872.d0 + l*(626.d0 + (-28.d0 + l)* l))))))))/xsq  &
               +  (22880.d0/256.d0) * (179928.d0 + l*(30408.d0 +  5.d0* l *(11596.d0 + l* (-429.d0 + l* (391.d0 + (-15.d0 + l)* l))))))/xsq &
               - (45760.d0/256.d0) * (1620.d0 + l* (54.d0 + l* (191.d0 + (-6.d0 + l)*l))))/xsq +  (7040.d0/256.d0)* (52.d0 + lsq-l))/x - x) *SphericalBesselJ(l-1,x) &
               + ((((((46189.d0/256.d0)*(1.d0 + l)*(2.d0 + l)*(3.d0 + l)*(4.d0 + l)*(5.d0 + l)*(6.d0 + l)*(7.d0 + l)*(8.d0 + l)*(9.d0 + l)* (10.d0 + l)/xsq -  (24310/256.d0)*(3265920.d0 +  l* (3360528.d0 + l* (3086460.d0 +  l*(957492.d0 +  5* l*(59233.d0 + l*(7272.d0 + l*(850.d0 + l*(36.d0 + l)))))))))/xsq +  (22880.d0/256.d0)* (211680.d0 +  l* (100212.d0 +  5.d0* l *(15394.d0 + l *(2163.d0 + l*(481.d0 + l* (21.d0 + l)))))))/xsq - (45760.d0/256.d0)*(1764.d0 + l* (350.d0 + l*(215.d0 + l*(10.d0 + l)))))/xsq + (7040.d0/256.d0)*(54.d0 + l*(3.d0 + l)))/x - x) *SphericalBesselJ(l+1,x) &
               )/(2.d0*l+1.d0)
       case default
          jllpm = jllp0(l, lp, x)
       end select
       return
    case(1) !!m=1
       select case(lp)
       case(1)
          jllpm = (sqrt_table(l)*sqrt_table(l+1)/const_sqrt2)*SphericalBesselJ(l,x)/x
       case(2)
          jllpm = ((const_sqrt3/const_sqrt2)*sqrt_table(l)*sqrt_table(l+1))*(&
               (l-1.d0)*SphericalBesselJ(l-1,x) &
               -(l+2.d0)*SphericalBesselJ(l+1,x) &
               )/x/(2.d0*l+1.d0)
       case(3)
          jllpm = (const_sqrt3*sqrt_table(l)*sqrt_table(l+1))*( &
               (1.25d0*(l-2)*(l-1.d0)/xsq-1.d0)*SphericalBesselJ(l-1,x) &
               +(1.25d0*(l+2.d0)*(l+3.d0)/xsq-1.d0)*SphericalBesselJ(l+1,x) &
               )/(2.d0*l+1.d0)
       case(4)
          jllpm = (const_sqrt5*sqrt_table(l)*sqrt_table(l+1))/x*( &
               (1.75d0*(l-3)*(l-2)*(l-1.d0)/xsq-(l-8.d0))*SphericalBesselJ(l-1,x) &
               - (1.75d0*(l+2.d0)*(l+3.d0)*(l+4.d0)/xsq-(l+9.d0))*SphericalBesselJ(l+1,x) &
               )/(2.d0*l+1.d0)
       case(5)
          jllpm = (sqrt(7.5d0)*sqrt_table(l)*sqrt_table(l+1))*( &
               (((21.d0/8.d0)*(l-4)*(l-3)*(l-2)*(l-1.d0)/xsq-3.5d0*(20.d0+lsq-3.d0*l))/xsq+1.d0)*SphericalBesselJ(l-1,x) &
               +(((21.d0/8.d0)*(l+2.d0)*(l+3.d0)*(l+4.d0)*(l+5.d0)/xsq-3.5d0*(24.d0+lsq+5.d0*l))/xsq+1.d0)*SphericalBesselJ(l+1,x) &
               )/(2.d0*l+1.d0)
       case(6)
          jllpm = (sqrt(10.5d0)*sqrt_table(l)*sqrt_table(l+1))/x*( &
               (((33.d0/8.d0)*(-5 + l)*(-4 + l)* (-3 + l)*(-2 + l)* (-1.d0 + l)/xsq-4.5d0*(-160.d0 + l*(22 + (-17.d0 + l)*l)))/xsq+l-19.d0)*SphericalBesselJ(l-1,x) &
               +(((-33.d0/8.d0)*(2.d0 + l)*(3.d0 + l)*(4.d0 + l)*(5.d0 + l)*(6.d0 + l)/xsq+4.5d0*(200.d0 + l*(59.d0 + l*(20.d0 + l))))/xsq-(l+20.d0))*SphericalBesselJ(l+1,x) &
               )/(2.d0*l+1.d0)
       case default
          jllpm = 0
       end select
       return
    case(2) !!m=2
       select case(lp)
       case(2)
          jllpm = (sqrt(0.375d0)*sqrt_table(l)*sqrt_table(l-1)*sqrt_table(l+1)*sqrt_table(l+2))*SphericalBesselJ(l,x)/xsq
       case(3)
          jllpm = (sqrt(15.d0/8.d0)*sqrt_table(l)*sqrt_table(l-1)*sqrt_table(l+1)*sqrt_table(l+2))/xsq*( &
               (l-2.d0)*SphericalBesselJ(l-1, x) &
               - (l+3.d0)*SphericalBesselJ(l+1, x) &
               )/(2.d0*l+1.d0)
       case(4)
          jllpm = (sqrt(2.5d0)*sqrt_table(l)*sqrt_table(l-1)*sqrt_table(l+1)*sqrt_table(l+2))/x*( &
               (1.75d0*(l-3)*(l-2)/xsq-1.5d0)*SphericalBesselJ(l-1,x) &
               +(1.75d0*(l+3.d0)*(l+4.d0)/xsq-1.5d0)*SphericalBesselJ(l+1,x) &
               )/(2.d0*l+1.d0)
       case(5)
          jllpm = (sqrt(105.d0/8.d0)*sqrt_table(l)*sqrt_table(l-1)*sqrt_table(l+1)*sqrt_table(l+2))/xsq*( &
               (1.5d0*(l-4)*(l-3)*(l-2)/xsq-(l-11.d0))*SphericalBesselJ(l-1,x) &
               + (-1.5d0*(l+3.d0)*(l+4.d0)*(l+5.d0)/xsq +(l+12.d0))*SphericalBesselJ(l+1,x) &
               )/(2.d0*l+1.d0)
       case(6)
          jllpm = (sqrt(105.d0/1024.d0)*sqrt_table(l)*sqrt_table(l-1)*sqrt_table(l+1)*sqrt_table(l+2))/x*( &
               ((33.d0*(l-5.d0)*(l-4.d0)*(l-3.d0)*(l-2.d0)/xsq-48.d0*(39.d0+l*(l-5.d0)))/xsq+16.d0)*SphericalBesselJ(l-1,x) &
               +((33.d0*(l+3.d0)*(l+4.d0)*(l+5.d0)*(l+6.d0)/xsq-48.d0*(45.d0+l*(l+7.d0)))/xsq+16.d0)*SphericalBesselJ(l+1,x) &
               )/(2.d0*l+1.d0)
       case default
          jllpm = 0
       end select
       return
    end select
  end function Jllpm

  function jllpm_vec(l,lp,m, x)
    real(dl),dimension(:),intent(in)::x
    real(dl) jllpm_vec(size(x))
    integer i, l, lp, m
    do i=1,size(x)
       jllpm_vec(i) = jllpm(l,lp,m,x(i))
    enddo
  end function jllpm_vec

  function jllp0_inner_product(l, lp, maxl, n, jl, y) result(ss)
    integer i, l, lp, maxl, n, jmin, j, lbase, ltop
    real(dl) jl(0:maxl, n), y(n), coef(0:lp)
    real(dl) ss
    ss = 0._dl
    jmin = max(0, lp-l)
    lbase = l-lp+jmin*2
    ltop = l+lp
    do j = jmin, lp
       coef(j)= jllp0_coef(l,lp,j)
    enddo
    !omp parallel do reduction(+:ss)
    do i=1, n
       if(abs(y(i)).gt. 1.e-15_dl)then
          ss = ss + sum(coef(jmin:lp)*jl(lbase:ltop:2, i)) *y(i)
       endif
    enddo
    !omp end parallel do
  end function jllp0_inner_product

  function jllp0_coef(l,lp,j) result(coef)
    real(dl) coef
    integer l, lp, j
    coef =(1.d0+2.d0*l+4.d0*j-2.d0*lp)*dexp( &
         -1.5d0*const_lnpi+(2*j-1)*const_ln2 &
         + lnGamma_int_table(l + j + 1) &
         + lnGamma_hi_table(1 + l + j-lp) &
         + lnGamma_hi_table(1 - j + lp) &
         - lnGamma_int_table(1 + 2*j) &
         - lnGamma_int_table(l + j - lp + 1) &
         - lnGamma_int_table(lp - j + 1) &
         + 2.d0 * lnGamma_hi_table(1 + j) &
         - lnGamma_hi_table(2 + l + j) &
         )
    if(mod(j,2).ne.0) coef = -coef
  end function jllp0_coef
!!$


  function jllp0(l, lp, x)
    real(dl) x, jllp0
    integer l, lp, j
    jllp0 = 0.d0
    do j=max(0, lp-l), lp
       jllp0 = jllp0 + jllp0_coef(l,lp,j)*SphericalBesselJ(l-lp+j*2,x)
    enddo
  end function jllp0


  subroutine get_jl_smoothw(lp, m, n, x, src, ns, w,  lmin, lmax)
    integer,intent(IN):: lp, m, ns, n, lmin, lmax
    real(dl),intent(IN):: x(n), src(n)
    real(dl),intent(OUT)::w(-ns:ns)
    real(dl) jl(lmin:lmax, n), jlm(lmin:lmax, ns+1:n-ns)
    real(dl) fx(lmin:lmax+1, -ns:ns), y(lmin:lmax+1)
    integer l, i
    !$omp parallel do private(i, l)
    do l = lmin, lmax
       do i=1, n
          jl(l, i) = SphericalBesselJ(l, x(i))
       enddo
       do i=ns+1, n-ns
          jlm(l, i) = jllpm(l, lp, m, x(i))
       enddo
    enddo
    !$omp end parallel do 
    do l=lmin, lmax
       do i=-ns,ns
          fx(l, i) = sum(jl(l, 1+ns-i:n-ns-i)*src(1+ns:n-ns)) 
       enddo
       if(m.gt.0) fx(l, -ns:ns)= fx(l,-ns:ns) * jllpm_lfact(l, m)
       y(l) = sum(jlm(l, 1+ns:n-ns)*src(1+ns:n-ns))
    enddo
    y(lmax+1) = 0.d0
    fx(lmax+1,-ns:ns) = 0.1d0**(lp-m)
    call fit_template(lmax-lmin+2, 2*ns+1, y(lmin:lmax+1), fx(lmin:lmax+1, -ns:ns), w(-ns:ns))
  end subroutine get_jl_smoothw


  function jllpm_lfact(l, m) result(fact)
    integer l, m, i
    real(dl) fact
    fact = 1.d0
    do i= l-m+1, l+m
       fact = fact * sqrt_table(i)
    end do
  end function jllpm_lfact


end module jllpm_utils

