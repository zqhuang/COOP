module smooth_kessence

contains

  function kessence_F(x) result(F)
    real*8::x, F, y
    y = x**3
    if(x > 0.03d0)then
       F = sqrt(1.d0+1.d0/y) - log(sqrt(y) + sqrt(1.d0+y))/y
    else
       F = ((2.d0/3.d0) + y*(- 0.2d0+y*(3.d0/28.d0)))*sqrt(y)
    endif
  end function kessence_F

  function kessence_F2(x) result(F2)  !!F^2(x)/x
    real*8 x, F2
    if(x > 2.d-4)then
       F2 = kessence_F(x) ** 2 / x
    else
       F2 = 4.d0/9.d0*x**3
    endif
  end function KESSENCE_F2


  function kessence_intF2(x) result(intF2)  !!\int_0^x  F^2(t)/t dt
    real*8,parameter::xpiv = 1.4d0
    real*8 x, intF2, y
    if(x < 0.801d0)then
       y = x**3       
       intF2 = y*(4.d0/27.d0 &
            + y * ( - 2.d0/45.d0 &
            + y * ( 32.d0/1575.d0 &
            + y * (-32.d0/2835.d0 &
            + y * ( 512.d0/72765.d0 &
            + y * (-128.d0/27027.d0 ))))))
    elseif(x < 3.d0) then
       y = x- xpiv
       intF2 =  0.24680697804518203    &
            + y * ( 0.35215257966981595   &
            + y * (        3.9796114258954955E-002  &
            + y * (   -6.2970564003774390E-002  &
            + y * (    2.6841690599631041E-002  &
            + y * (    8.1426587339724753E-004  &
            + y * (   -8.1468541931120960E-003  &
            + y * (    4.3342792297875413E-003  &
            + y * (   -7.9148050497320540E-004))))))))
    else
       intF2 = log(x) -0.37120535360415519  &
            + (-2.5375113621010333E-003 &
            + (0.10438085428501427      &
            + (  2.3611450907305018     &
            + (  -4.5758725635158362    &
            + (    3.8113837616994544   &
            + (    0.65371487085641022  &
            + (     -2.1964577462867707 &
            + (-2.4105961820664401/x))/x)/x)/x)/x)/x)/x)/x
    endif
  end function KESSENCE_INTF2
  
  function kessence_dFdx(x) result(dFdx)
    real*8::x, dFdx, y,sqrty, sqrtz
    y = x**3
    if(x > 0.03d0)then
       sqrty = sqrt(y)
       sqrtz = sqrt(1.d0+y)
       dFdx = 3.d0/x/sqrty*(log(sqrty+sqrtz)/sqrty-1.d0/sqrtz)
    else
       dFdx = (1.d0 + y*( - 0.9d0 + y* (45.d0/56.d0)))*sqrt(x)
    endif
  end function kessence_dFdx


  function wrat(x)  !! w_a/(1+w_0) as a function of 1/a_eq
    real*8 x, wrat, y
    if(x > 0.3d0)then
       wrat = -2*x*kessence_dFdx(x)/kessence_F(x)
    else
       y = x**3
       wrat = -3.d0+y*(1.8d0+y*(-243.d0/175.d0 + y * (1023.d0/875.d0- y*(346431.d0/336875.d0))))
    end if
  end function wrat


  function inverse_wrat(wr) result(x)
    real*8 wr, x, c, diff
    real*8,parameter::a = 243.d0/175.d0, b=1.8d0
    real*8,parameter::cc = 6.d0*(1.d0-log(2.d0)), c2= - 15.d0/2 - 6*log(2.d0)**2 + 9*log(2.d0)
    real*8,parameter::wrpiv = -1.3d0
    c = wr + 3.d0
    if(c < 1.d-5)then
       x = (c/b)**(1.d0/3.d0)
    else if(c < 0.21 ) then
       x = ((b - sqrt(b**2-4*a*c))/(2*a))**(1.d0/3.d0)
    else if( wr < - 0.33 )then
       c = wr - wrpiv
       x = 1.4144315993670178  &
            + c*(  0.77052678129513374   &
            + c*( 0.26126811863262378    &
            + c*( 0.11198575034585553    &
            + c*( 0.19960058705922162    &
            + c*( 0.25444465810140654    &
            + c*( 8.4157349326106523E-002))))))
    else
       x = 5.d0
       diff =  wr - asymp(x)
       do while( abs(diff)> abs(wr)*1.e-4)
          x = x + diff /  asymp_derv(x)
          diff =  wr - asymp(x)          
       end do
    endif
  contains
    function asymp(x)
      real*8 x, asymp
      asymp = (cc-9*log(x)+ ((18*log(2.d0)-13.5d0-13.5d0*log(x))*log(x) + c2)/x**3)/x**3
    end function asymp

    function asymp_derv(x)
      real*8 x, asymp_derv
      asymp_derv = -3.d0*(cc-9.d0*log(x)-9.d0/x)/x**4
    end function asymp_derv
  end function inverse_wrat


  function kessence_wofa(a, w0, wa) result(w)
    real*8 a, w0, wa , w, invaeq
    if(1.d0+w0 .eq. 0.d0)then
       w = -1.d0
    else
       invaeq = inverse_wrat(wa/(1.d0+w0))
       w = -1.d0 + (1.d0+w0) * (kessence_F(a*invaeq)/kessence_F(invaeq))**2
    endif
  end function kessence_wofa

  function kessence_rhoofa(a, w0, wa) result(rho) !! return rho(a)/rho(z=0)
    real*8 a, w0, wa , rho, invaeq
    if(1.d0+w0 .eq. 0.d0)then
       rho = 1.0
    else
       invaeq = inverse_wrat(wa/(1.d0+w0))
       rho = exp(3.d0*(1.d0+w0)/kessence_F(invaeq)**2*(kessence_intF2(invaeq)-kessence_intF2(a*invaeq)))
    endif
  end function kessence_rhoofa


  function kessence_wa(w0, Omega_m)result(wa)
    real*8::w0, Omega_m, wa
    wa = (1.d0+w0)*wrat(((1.d0-Omega_m)/Omega_m)**(1.d0/3.d0))
  end function kessence_wa

end module smooth_kessence



program Test
  use smooth_kessence
  use coop_wrapper_utils
  implicit none
  integer,parameter::n=128
  real*8:: Omega_m(n), wr(n), omm
  integer i
  type(coop_asy)::fig
  print*,coop_LegendreP(10000, cos(coop_pi/100.))
  print*, sqrt(2/coop_pi/(100*coop_pi))/sqrt(2.)
  stop
  call fig%open("wrat.txt")
  call fig%init(xlabel = "$\Omega_m$", ylabel = "$\frac{w_a}{1+w_0}$")
  do i = 1, n
     Omega_m(i) = 0.01+0.98*(i-1)/(n-1.d0)
     wr(i) =  wrat(((1.d0-Omega_m(i))/Omega_m(i))**(1.d0/3.d0))
  enddo
  call fig%plot(Omega_m, wr, color="blue", linewidth=2., legend="$\frac{w_a}{1+w_0}=\frac{-2\left(\frac{1-\Omega_m}{\Omega_m}\right)^{1/3}F'\left[\left(\frac{1-\Omega_m}{\Omega_m}\right)^{1/3}\right]}{F\left[\left(\frac{1-\Omega_m}{\Omega_m}\right)^{1/3}\right]}$")
  wr = -1.42*(Omega_m/0.3)**0.64
  call fig%plot(Omega_m, wr, color="orange", linetype = "dotted", linewidth=2., legend="$\frac{w_a}{1+w_0}=  -1.42 \left(\frac{\Omega_m}{0.3}\right)^{0.64}$")
  call fig%legend(0.4, 0.9)
  call fig%close()
end program Test



