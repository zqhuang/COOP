module coop_sbj_mod
  use coop_wrapper_typedef
  implicit none

#include "constants.h"

  private
  integer,parameter::jlb_lmax = 5101
  integer,parameter::jlb_lmin = 10
  integer,parameter::dl = kind(1.d0)

  type jlb_saved_data
     real(dl) dx, xmax, xmin
     integer n
     real(dl),dimension(:,:),allocatable::data
  end type jlb_saved_data

  logical::jlb_tabulated(jlb_lmin:jlb_lmax) = .false.
  type(jlb_saved_data)::jlb_table(jlb_lmin:jlb_lmax)

!!yes -> 1.3 GB mem
!!no ->  0.4 GB mem
#define  JLB_SUPER_PRECISION COOP_NO 

contains


  function jlb_phase(l, x) result(phase)
    integer l
    real(dl) phase
    real(dl) x,  nu2, nu, chb, ch2b, sh2b, shb, sx, sat, rat
    if(l.eq.0)then
       phase  = x - coop_pio2
       return
    endif
    nu2 = l*(l+1.d0)
    nu = sqrt(nu2)
    if(x .le. nu)then
       phase = 0
       return
    endif
    chb = x/nu
    ch2b = x**2/nu2
    sh2b = ch2b - 1.d0
    shb = sqrt(sh2b)
    sx = shb * nu
    sat = 1.d0/shb/sh2b/nu
    rat = sat*sat
    phase = (-coop_pio4 &
         + sx &
         -((-1.d0/128/nu2 + 0.125)/nu +nu)*acos(nu/x) &
         + sat*( &
         -1.d0/12.d0-ch2b/8.d0  &
         + rat * (-31/5760.d0 + ch2b *( 1062/5760.d0 + ch2b * (4644/5760.d0+ ch2b * ( -195/5760.d0 + 45/5760.d0 * ch2b))) &
         - rat * ch2b**2* (174126/32256.d0 + ch2b*(511445/32256.d0+ ch2b*(180144/32256.d0)) &
         ))))
  end function jlb_phase

  function sphericalbesselj_zero(l, k) result(z)
    integer k, l
    real(dl) z, t, zmin, zmax, phasemax, phase, phasemid
    real(dl),parameter::c1(12) =   (/ -0.245127043d+01, 0.112632772d+02, -0.207925322d+02, 0.198543806d+02, -0.110527677d+02, 0.406769876d+01, -0.109067653d+01, 0.206080657d+00, -0.260196064d-1 , 0.212880935d-2 , -0.101897908d-3 , 0.216576354d-5  /)
    real(dl), parameter::c2(12) =  (/ 0.579779690d+01, -0.280853646d+02, 0.625315957d+02, -0.805135267d+02, 0.629961403d+02, -0.307235514d+02, 0.973806423d+01, -0.206907154d+01, 0.293355044d+00, -0.265358216d-1 , 0.138530622d-2 , -0.317522214d-4  /)
    if(l.eq.0)then
       z = coop_pi*k
       return
    endif
    select case(k)
    case(1)
       t = (l+0.5d0)**0.25d0
       z = l + 0.5d0 + 2.d0*(l+2.5d0)**0.323d0  + polyvalue(c1, t)/(l+0.5d0)
    case(2)
       t = (l+0.5d0)**0.25d0
       z = l + 0.5d0 + 3.51d0*(l+5.d0)**0.323d0 + polyvalue(c2, t)/(l+0.5d0)
    case default
       t = (l+0.5d0)**0.25d0
       zmin = l + 0.5d0 + 3.51d0*(l+5.d0)**0.323d0 + polyvalue(c2, t)/(l+0.5d0)
       phase = coop_pi * (k-0.5_dl) 
       zmax = max(sphericalbesselj_zero_asym(l, k), zmin+coop_pi)
       phasemax = jlb_phase(l, zmax)
       do while(phasemax .lt. phase)
          zmax = zmax + max(phase - phasemax, coop_pi)
          phasemax = jlb_phase(l, zmax)
       enddo
       do while(zmax - zmin .gt. 1.d-6)
          z = (zmax+zmin)/2.d0
          phasemid = jlb_phase(l, z)
          if(phasemid .ge. phase)then
             zmax = z
          else
             zmin = z
          endif
       enddo
       z = (zmin+zmax)/2.d0
    end select
  end function sphericalbesselj_zero


!!the k-th zero of spherical Bessel function j_l
  function sphericalbesselj_zero_asym(l, k) result(z)
    real(dl) z, nusq, zsq
    integer l, k
    z = coop_pi * (l/2.d0+k)
    nusq = (l +0.5d0)**2
    zsq = z**2
    z = z - l*(l+1.d0)/(2.d0*z)*(1.d0+ ((-31.d0+28.d0*nusq) + &
         (3779.d0+ nusq*(-3928.d0 + 1328.d0*nusq)+(-6277237.d0+ nusq*(6342972.d0 + nusq*(-2461680.d0+nusq*444736.d0)))/(224.d0*zsq))/(40.d0*zsq))/(48.d0*zsq))
  end function sphericalbesselj_zero_asym



!!start where j_l is ~ 10^{-8}
  subroutine jlb_startpoint(l, xstart, eps)
    real(dl) xstart, lnnu, nu
    real(dl),optional::eps
    integer l
    nu = l+ 0.5d0
    lnnu = log(nu)
    select case(l)
    case(0:5)
       xstart = exp(-26.8614 + lnnu * ( 28.1288 - lnnu * 8.14055))
    case(6:20)
       xstart = nu - exp(-0.233555 + lnnu*(1.369 - lnnu* 0.145824))
    case default
       xstart = nu - exp(1.52077 -2.28973/nu+ lnnu*(0.41562 -0.00916079 * lnnu))
    end select
    if(present(eps) .and. l .gt. 0)then
       xstart = xstart*(eps/1.e-8_dl)**(1.d0/l)
    endif
  end subroutine jlb_startpoint



!!very accurate for x < l - (a few) * l^{1/3}
 subroutine jlb_small(l, x, jl, jlp)
    integer l
    real(dl),intent(IN)::x
    real(dl)nu, nu2, x2, cosb, cos2b, sinb, sin2b, sin3b,  nusin3b, invnu3b, expterm
    real(dl),intent(OUT):: jl
    real(dl),intent(OUT),optional:: jlp
    if(l.lt.jlb_lmin)then
       if(present(jlp))then
          call jlb_exact_lowl(l, x, jl, jlp)
       else
          call jlb_exact_lowl(l, x, jl)
       endif
       return
    endif
    x2 = x**2
    if(x2 .lt. l*3)then
       if(x .lt. 1.d-20) then
          jl = 0.d0
          if(present(jlp)) jlp = 0.d0
          return
       endif
       cosb = x2/4.d0
       expterm = dexp(-lngamma_hi_table(l+2)+l*(dlog(x)-coop_ln2)-coop_ln2)*coop_sqrtpi 
       jl =  expterm * (1.d0  - cosb /(l+1.5d0)*(1.d0 - cosb/(l+2.5d0)/2.d0*(1.d0 - cosb/(l+3.5d0)/3.d0*(1.d0-cosb/(l+4.5d0)/4.d0*(1.d0 - cosb/(l+5.5d0)/5.d0 * (1.d0 - cosb/(l+6.5d0)/6.d0*(1.d0-cosb/(l+7.5d0)/7.d0*(1.d0-cosb/(l+8.5d0)/8.d0*(1.d0-cosb/(l+9.5d0)/9.d0*(1.d0-cosb/(l+10.5d0)/10.d0*(1.d0-cosb/(l+11.5d0)/11.d0*(1.d0-cosb/(l+12.5d0)/12.d0))))))))))))
       if(present(jlp))then
          jlp = l/x * jl + expterm * (-x/ 4.d0 /(l+1.5d0)*(2.d0 - cosb/(l+2.5d0)/2.d0*(4.d0 - cosb/(l+3.5d0)/3.d0*(6.d0-cosb/(l+4.5d0)/4.d0*(8.d0 - cosb/(l+5.5d0)/5.d0 * (10.d0 - cosb/(l+6.5d0)/6.d0*(12.d0-cosb/(l+7.5d0)/7.d0*(14.d0-cosb/(l+8.5d0)/8.d0*(16.d0-cosb/(l+9.5d0)/9.d0*(18.d0-cosb/(l+10.5d0)/10.d0*(20.d0-cosb/(l+11.5d0)/11.d0*(22.d0-cosb/(l+12.5d0)*2.d0))))))))))))
       endif
       return
    endif
    nu2 =l*(l+1.d0)  
    nu = sqrt(nu2)
    x2 = x**2
    cosb = x/nu
    cos2b = cosb**2
    sin2b = 1.d0 - cos2b
    sinb = sqrt(sin2b)
    sin3b = sin2b * sinb
    nusin3b = nu * sin3b
    invnu3b = 1.d0/nusin3b
    jl = dexp( &
         nu*sinb - coop_ln2 &
         + invnu3b * ( -1.d0/12.d0 - 1.d0/8.d0 * cos2b &
         + invnu3b * ( -1.d0/16.d0 + 3.d0/8.d0 *cos2b &
         + invnu3b * (31.d0/5760 + cos2b*(-1062.d0/5760 + cos2b*( -4644.d0/5760 + cos2b*(195.d0/5760+ cos2b*(-45.d0/5760))))  &
         + invnu3b * (1.d0/128+cos2b*(12.d0/128+cos2b*(312.d0/128+ cos2b*(240.d0/128)))  &
         + invnu3b * ((-(15.9- 57 *invnu3b)*cos2b - 5.4)*cos2b**2 &
         ))))) &
         -0.25d0 * dlog(x2*(nu2-x2)) &
         + (l+0.5d0)*dlog(cosb/(1+sinb)) &
         )
    if(present(jlp))then
       jlp = jl * ( &
            x/nu*( &
            (( -3.d0/8.d0 * (1-6*cos2b)/sin2b + 0.75d0)/(nu*sin3b) &
            - 0.25d0 &
            - (0.25 +3.d0/8.d0*cos2b)/sin2b ) /sin3b/nu2 &
            -1.d0/sinb) &
            + ( &
            (24 + (1248  + 1440*cos2b)* cos2b)*cos2b/(128*x) &
            +(3*x*sinb * (1 +  cos2b*(12 + cos2b*(312 + cos2b*240 ))))/(32*nu * nusin3b) &
            + cos2b*(-2124+cos2b*(-18576 + cos2b*( 1170 - 360*cos2b )))*nusin3b/5760/x &
            + cosb * sinb*(31/640.d0 + cos2b*(- 1062/640.d0 +cos2b*(- 4644/640.d0 + cos2b*(195/640.d0 - cos2b*(45/640.d0)))))  &
            -(21.6 + cos2b*(95.4 -342/(nu*sin3b)+ cos2b * (- 171/(nu*sin2b*sin3b))))*cos2b**2/x/nusin3b &
            -(15*x*cos2b**2*(5.4 + cos2b* (15.9-57/nusin3b)))/(nu2*nusin3b*sin2b) &
            )/nusin3b**4 &
            +( &
            +(cos2b-0.5d0)/sin2b &
            +(((1.d0/1024/nu2) - 1.d0/128)/nu2 + 0.125d0 + nu2)*( (cos2b/(sinb*(1+sinb))+1.d0)/nu) &
            )/x &
            )
    endif
  end subroutine jlb_small



  subroutine jlb_large(l, x, jl, jlp)
    integer l
    real(dl),intent(IN)::x
    real(dl),intent(OUT)::jl
    real(dl),intent(OUT),optional::jlp
    real(dl) sx, expterm, alpha, x2, nu, nu2, sh2b, shb, ch2b, chb, sat, rat, expderv, alphaderv, cschb, csch2b, invnu2, invnu
    if(l.lt.jlb_lmin)then
       call jlb_exact_lowl(l, x, jl, jlp)
    else
       x2 = x*x
       nu2 = l*(l+1.d0)
       nu = sqrt(nu2)
       chb = x/nu
       ch2b = x2/nu2
       sh2b = ch2b - 1.d0
       shb = sqrt(sh2b)
       sx = shb * nu
       sat = 1.d0/shb/sh2b/nu
       rat = sat*sat
       expterm = dexp( &
            (1.d0/16-3.d0/8.d0*ch2b  &
            + ( 1.d0/128.d0 + ch2b*(3.d0/32.d0 + ch2b * (39/16.d0 + ch2b * (15.d0/8.d0))) &
            - (2523.d0/32.d0 + ch2b*(1521.d0/16.d0 + ch2b * (315.d0/16.d0)))*ch2b**3 * rat)*rat)*rat &
            -0.5d0*dlog(x*sx) &
            )
       alpha = -coop_pio4 &
            + sx &
            -((-1.d0/128/nu2 + 0.125)/nu +nu)*acos(nu/x) &
            + sat*( &
           -1.d0/12.d0-ch2b/8.d0  &
            + rat * (-31/5760.d0 + ch2b *( 1062/5760.d0 + ch2b * (4644/5760.d0+ ch2b * ( -195/5760.d0 + 45/5760.d0 * ch2b))) &
            - rat * ch2b**2* (174126/32256.d0 + ch2b*(511445/32256.d0+ ch2b*(180144/32256.d0)) &
            )))
       jl = expterm*dcos(alpha)
       if(present(jlp))then
          csch2b = 1.d0/sh2b
          cschb = 1.d0/shb
          invnu = 1.d0/nu
          invnu2 = 1.d0/nu2
          expderv = invnu*chb*csch2b*( &
               -0.5d0 + 3.d0/32.d0*invnu2*csch2b**2*( &
               16 &
               + csch2b * (20-120*invnu2 &
               + csch2b * (16*invnu2*(-43+105*invnu2) &
               + csch2b * (10*invnu2*(-113+2064*invnu2) &
               + csch2b * (invnu2*(-565+83964*invnu2) &
               + csch2b * (149898*invnu2**2 &
               + csch2b * (122064*invnu2**2 &
               + csch2b * (37170*invnu2**2 &
               ))))))))) &
               - 0.5d0/x


          alphaderv = (1/10752.d0)*( &
               invnu2*cschb*( &
               chb*csch2b*( &
               -84*(-16+invnu2) &
               + csch2b * ( 84*(80+invnu2) &
               + csch2b * (-40404*invnu2 &
               + csch2b * ( 84*invnu2*(-1547+5004*invnu2) &
               + csch2b * ( 21*invnu2*(-4420+176003*invnu2) &
               + csch2b * ( 10227525*invnu2**2 &
               + csch2b * ( 11280373*invnu2**2 &
               + csch2b * ( 4328575*invnu2**2 &
               )))))))) &
               +84*(-16+invnu2)/chb &
               ) ) &
               + shb/chb

          jlp = jl * expderv - expterm * dsin(alpha) * alphaderv
       endif
    end if
  end subroutine jlb_large

  subroutine jlb_extra_small(l, x, jl)
    integer l
    real(dl) x, jl, nu, nu2, sx, x2
    x2  = x*x
    if(x2.le.l*2)then
       if(x.lt. 1.d-30)then
          jl = 0.d0
          return
       endif
       x2 = x2/4.d0
       jl = dexp(-lngamma_hi_table(l+2)+l*(dlog(x)-coop_ln2)-coop_ln2)*coop_sqrtpi  * (1.d0  - x2 /(l+1.5d0)*(1.d0 - x2/(l+2.5d0)/2.d0*(1.d0 - x2/(l+3.5d0)/3.d0*(1.d0-x2/(l+4.5d0)/4.d0*(1.d0 - x2/(l+5.5d0)/5.d0 * (1.d0 - x2/(l+6.5d0)/6.d0*(1.d0-x2/(l+7.5d0)/7.d0*(1.d0-x2/(l+8.5d0)/8.d0))))))))
       return
    endif
    nu = l+0.5d0 - 0.125d0/(l+0.5d0)
    nu2 = nu*nu    
    sx = dsqrt(nu2-x2)
    jl =  dexp(-(l+0.5d0)*dlog((nu +sx)/x) + sx  - (1.d0/12.d0*nu2+1.d0/8.d0*x2)/sx**3 - coop_ln2 - dlog(sx*x)/2.d0)
  end subroutine jlb_extra_small

  subroutine jlb_extra_large(l, x, jl)
    integer l
    real(dl) x, jl
    real(dl) sx, nu, nu2, rat, sat, nubyx
    nu = l+0.5d0 
    nu2 = nu*nu
    nubyx = nu/x
    rat = nubyx**2
    sat = 1.d0 - rat * ( 0.25d0 + rat * ( 3.d0/32.d0 + rat * ( 7.d0/128.d0 + rat * ( 77.d0/2048.d0 + rat * ( 231.d0/8192.d0 + rat * ( 1463.d0/65536.d0 + rat * ( 4807.d0/262144.d0 + rat * (129789.d0/8388608.d0 + rat * (447051.d0/33554432.d0 + rat*(3129357.d0/268435456.d0)))))))))) 
    sx = x * sat**2
    jl = ( 1.d0 - ( 0.25d0*nu2 + x**2/16. ) / sat**2/sx**4 ) * dcos( &
         sx - nu * ( coop_pio2 - nubyx * ( 1.d0 + rat * (1.d0/6.d0 + rat * (3.d0/40.d0 + rat * (5.d0/112.d0 + rat * (35.d0/1152.d0 + rat * (63.d0/2816.d0 + rat * (231.d0/13312.d0 + rat * (143.d0/10240.d0+rat*(6435.d0/557056.d0) ) ) ) ) ) ) ) ) ) - coop_pio4 &
         - ((1.d0/8.d0)/rat+1.d0/12.d0)*nu2/sx**3 &
         ) / x / sat
  end subroutine jlb_extra_large

  subroutine jlb_exact_lowl(l, x, jl, jlp)
    integer l
    real(dl) x, jl, v2, x2
    real(dl),optional::jlp
    if(x .lt. l/7.d0+0.1d0)then
       x2 = x**2
       select case(l)
       case(0)
          jl = 1.d0 + x2*(-1.d0/6.d0 + x2*(1.d0/120.d0 + x2*(-1.d0/5040.d0 + x2*(1.d0/362880.d0+x2*(-1.d0/39916800.d0)))))
          if(present(jlp))jlp=x*(-1.d0/3.d0+x2*(1.d0/30.d0+x2*(-1.d0/840.d0+x2*(1.d0/45360.d0+x2*(-1.d0/3991680.d0)))))
       case(1)
          jl = x*(1.d0/3.d0 + x2*(-1.d0/30.d0 + x2*(1.d0/840.d0 + x2*(-1.d0/45360.d0 + x2 * (1.d0/3991680.d0)))))
          if(present(jlp))jlp=1.d0/3.d0 + x2*(-0.1d0 + x2*(168.d0 + x2*(-1.d0/6480.d0)))
       case(2)
          jl = x2*(1.d0/15.d0 + x2*(-1.d0/210.d0 + x2*(1.d0/7560.d0 + x2*(-1.d0/498960.d0 + x2*(1.d0/51891840.d0)))))
          if(present(jlp))jlp=x*(2.d0/15.d0 + x2*(-2.d0/105.d0+x2*(1.d0/1260.d0+x2*(-1.d0/62370.d0+x2*(1.d0/5189184.d0)))))
       case(3)
          jl = x*x2*(1.d0/105.d0 + x2*(-1.d0/1890.d0+ x2*(1.d0/83160.d0 + x2*(-1.d0/6486480.d0))))
          if(present(jlp))jlp= x2*(1.d0/35.d0+x2*(-1.d0/378.d0+x2*(1.d0/11880.d0+x2*(-1.d0/720720.d0))))
       case(4)
          jl = x2*x2*(1.d0/945.d0 + x2*(-1.d0/20790.d0 + x2 * (1.d0/1081080.d0 + x2 *(-1.d0/97297200.d0))))
          if(present(jlp))jlp=x2*x*(4.d0/945.d0 + x2*(-6.d0/20790.d0 + x2 * (8.d0/1081080.d0 + x2 *(-1.d0/9729720.d0))))
       case(5)
          jl = x2*x2*x*(1.d0/10395.d0 + x2*(-1.d0/270270.d0 + x2*(1.d0/16216200.d0+x2*(-1.d0/1654052400.d0))))
          if(present(jlp))jlp=x2*x2*(5.d0/10395.d0 + x2*(-7.d0/270270.d0 + x2*(9.d0/16216200.d0+x2*(-11.d0/1654052400.d0))))
       case(6)
          jl = x2**3*(1.d0/135135.d0 + x2*(-1.d0/4054050.d0 + x2 * (1.d0/275675400.d0+x2*(-1.d0/31426995600.d0))))
          if(present(jlp))jlp= x2**2*x*(6.d0/135135.d0 + x2*(-8.d0/4054050.d0 + x2 * (1.d0/27567540.d0+x2*(-12.d0/31426995600.d0))))
       case(7)
          jl = x2**3*x*(1.d0/2027025.d0+x2*(-1.d0/68918850.d0+x2*(1.d0/5237832600.d0)))
          if(present(jlp))jlp=x2**3*(7.d0/2027025.d0+x2*(-9.d0/68918850.d0+x2*(11.d0/5237832600.d0)))
       case(8)
          jl = x2**4*(1.d0/34459425.d0+x2*(-1.d0/1309458150.d0+x2*(1.d0/109994484600.d0)))
          if(present(jlp))jlp=x2**3*x*(8.d0/34459425.d0+x2*(-1.d0/130945815.d0+x2*(12.d0/109994484600.d0)))
       case(9)
          jl = x2**4*x*(1.d0/654729075.d0+x2*(-1.d0/27498621150.d0+x2*(1.d0/2529873145800.d0)))
          if(present(jlp))jlp= x2**4*x*(9.d0/654729075.d0+x2*(-11.d0/27498621150.d0+x2*(13.d0/2529873145800.d0)))
       case default
          print*, "l = ", l
          stop "jlb_exact_lowl called for high l"
       end select
       return
    endif
    select case(l)
    case(0)
       jl = dsin(x)/x
       if(present(jlp))jlp = (dcos(x)-jl)/x
    case(1)
       jl = (dsin(x)/ x - dcos(x))/x
       if(present(jlp))jlp = (dsin(x) - 2.d0*jl)/x
    case(2)
       jl = ( dsin(x) * (-1.d0+ 3.d0/x**2) + dcos(x) * (-3.d0)/ x )/x
       if(present(jlp))jlp = (dsin(x)/ x - dcos(x) - 3.d0*jl)/x
    case(3)
       x2 = x*x
       jl = ( dsin(x) * (-6.d0+ 15.d0/x2)/ x & 
            + dcos(x) * ( 1.d0 -15.d0/x2))/x
       if(present(jlp))jlp = (dsin(x) * (-1.d0+ 3.d0/x**2) + dcos(x) * (-3.d0)/ x - 4.d0*jl)/x
    case(4)
       v2 = 1.d0/x**2
       jl = ( dsin(x) * ( &
            1.d0+ v2*(-45.d0+ v2*(105.d0))) &
            + dcos(x) * ( &
            10.d0 + v2*(-105.d0))/ x &
            )/x
       if(present(jlp))jlp = ( dsin(x) * (-6.d0+ 15.d0*v2)/ x & 
            + dcos(x) * ( 1.d0 -15.d0*v2) - 5.d0*jl)/x
    case(5)
       v2 = 1.d0/(x*x)
       jl = ( dsin(x) * ( &
            15.d0+ v2*(-420.d0+ v2*(945.d0)))/ x &
            + dcos(x) * ( &
            -1.d0 + v2*(105.d0 + v2*(-945.d0))) &
            )/x
       if(present(jlp))jlp = ( dsin(x) * ( &
            1.d0+ v2*(-45.d0+ v2*(105.d0))) &
            + dcos(x) * ( &
            10.d0 + v2*(-105.d0))/ x - 6.d0*jl)/x
    case(6)
       v2 = 1.d0/(x*x)
       jl = ( dsin(x) * ( &
            -1.d0+ v2*(210.d0+ v2*(-4725.d0+ v2*(10395.d0)))) &
            + dcos(x) * ( &
            -21.d0 + v2*(1260.d0 + v2*(-10395.d0)))/ x &
            )/x
       if(present(jlp))jlp = (  dsin(x) * ( &
            15.d0+ v2*(-420.d0+ v2*(945.d0)))/ x &
            + dcos(x) * ( &
            -1.d0 + v2*(105.d0 + v2*(-945.d0))) &
            -7.d0*jl)/x
    case(7)
       v2 = 1.d0/(x*x)
       jl = ( dsin(x) * ( &
            -28.d0+ v2*(3150.d0+ v2*(-62370.d0+ v2*(135135.d0))))/ x &
            + dcos(x) * ( &
            1.d0 + v2*(-378.d0 + v2*(17325.d0 + v2*(-135135.d0)))) &
            )/x
       if(present(jlp)) jlp = ( dsin(x) * ( &
            -1.d0+ v2*(210.d0+ v2*(-4725.d0+ v2*(10395.d0)))) &
            + dcos(x) * ( &
            -21.d0 + v2*(1260.d0 + v2*(-10395.d0)))/ x &
            -8.d0*jl)/x 
    case(8)
       v2 = 1.d0/(x*x)
       jl = ( dsin(x) * ( &
            1.d0+ v2*(-630.d0+ v2*(51975.d0+ v2*(-945945.d0+ v2*(2027025.d0))))) &
            + dcos(x) * ( &
            36.d0 + v2*(-6930.d0 + v2*(270270.d0 + v2*(-2027025.d0))))/ x &
            )/x
       if(present(jlp)) jlp = (dsin(x) * ( &
            -28.d0+ v2*(3150.d0+ v2*(-62370.d0+ v2*(135135.d0))))/ x &
            + dcos(x) * ( &
            1.d0 + v2*(-378.d0 + v2*(17325.d0 + v2*(-135135.d0)))) &
            -9.d0*jl)/x
    case(9)
       v2 = 1.d0/(x*x)
       jl = ( dsin(x) * ( &
            45.d0+ v2*(-13860.d0+ v2*(945945.d0+ v2*(-16216200.d0+ v2*(34459425.d0)))))/ x &
            + dcos(x) * ( &
            -1.d0 + v2*(990.d0 + v2*(-135135.d0 + v2*(4729725.d0 + v2*(-34459425.d0))))) &
            )/x
       if(present(jlp)) jlp = ( dsin(x) * ( &
            1.d0+ v2*(-630.d0+ v2*(51975.d0+ v2*(-945945.d0+ v2*(2027025.d0))))) &
            + dcos(x) * ( &
            36.d0 + v2*(-6930.d0 + v2*(270270.d0 + v2*(-2027025.d0))))/ x &
            -10.d0*jl)/x
    case default
          print*, "l = ", l
          stop "jlb_exact_lowl called for high l"
    end select
  end subroutine jlb_exact_lowl


  subroutine jlb_tabulate(l, max_x)
    real(dl),optional::max_x
    integer l, nsteps, ii , jj, i, j
    real(dl) xmax, step, y(3), lsq, scale1, scale2, rl
    integer, parameter :: s = 3
    real(dl)::step_max 
    integer recur_steps 
    real(dl), parameter ::   b(s) = (/ 5.d0/18.d0, 4.d0/9.d0, 5.d0/18.d0/)
    ! Butcher tableau for 6th order Gauss-Legendre method
    real(dl), parameter :: a(s,s) =  reshape( (/ &
         5.d0/36.d0, 2.d0/9.d0 - 1.d0/sqrt(15.d0), 5.d0/36.d0 - 0.5d0/sqrt(15.d0), &
         5.d0/36.d0 + sqrt(15.d0)/24.d0, 2.d0/9.d0, 5.d0/36.d0 - sqrt(15.d0)/24.d0, &
         5.d0/36.d0 + 0.5d0/sqrt(15.d0), 2.d0/9.d0 + 1.d0/sqrt(15.d0), 5.d0/36.d0 /), (/ s, s /) )
    real(dl) g(3,s)
    if(l.lt.jlb_lmin)return
    if(present(max_x))then
       xmax = max_x
    else
#if JLB_SUPER_PRECISION
       xmax = min(max(l*2.8d0, 600.d0), l+2500.d0)
#else
       xmax = min(max(l*2.2d0, 300.d0), l+1000.d0)   
#endif
    endif
    if(xmax .le. 0 .or.  l.gt. jlb_lmax)then
       write(*,*) "jlb_tabulate: wrong input"
       stop
    endif
    if(jlb_tabulated(l))then
       if(xmax .le. jlb_table(l)%xmax)return
       deallocate(jlb_table(l)%data)
    endif
#if JLB_SUPER_PRECISION
    step_max = 0.02d0 + l/12900.d0
#else
    step_max = 0.031d0 + (l/6800.d0)
#endif
    recur_steps = 7 + floor(step_max/0.08)
    jlb_tabulated(l) = .true.
    jlb_table(l)%xmax = xmax
    call jlb_startpoint(l, jlb_table(l)%xmin, 1.d-24)
    jlb_table(l)%n = max(min(16384, ceiling((jlb_table(l)%xmax-jlb_table(l)%xmin)/(step_max))),64)
    jlb_table(l)%dx = (jlb_table(l)%xmax-jlb_table(l)%xmin)/jlb_table(l)%n
    scale1 = jlb_table(l)%dx
    scale2 = scale1**2
    allocate(jlb_table(l)%data(0:2, 0:jlb_table(l)%n))
    lsq = l*(l+1.d0)
    y(1) = jlb_table(l)%xmin
    call jlb_small(l, y(1), y(2), y(3))
    call fcn(y, g(:,1))
    ii = 0
    jlb_table(l)%data(0,ii) = y(2)
    jlb_table(l)%data(1,ii) = y(3)*scale1 
    jlb_table(l)%data(2,ii) = g(3,1)*scale2

    nsteps = ceiling(jlb_table(l)%dx/ min((jlb_table(l)%xmin+ii*jlb_table(l)%dx)/l, 1.d0)/min(step_max, 0.15d0))
    step = jlb_table(l)%dx/nsteps
    rl = l*1.3d0
    do while(ii .lt. jlb_table(l)%n)
       call do_step()
       if(y(1) .ge. rl )exit
       nsteps = ceiling(jlb_table(l)%dx/(y(1)/rl)/min(step_max,0.15d0))
       step = jlb_table(l)%dx/nsteps
    enddo

    do while(ii .lt. jlb_table(l)%n)
       call do_step()
    enddo

  contains

    subroutine fcn(yy, yyp)
      real(dl) yy(3), yyp(3)
      yyp(1) = 1.d0
      yyp(2) = yy(3)
      yyp(3) =  - ((yy(1) - lsq/yy(1))*yy(2) + 2 * yy(3))/yy(1)
    end subroutine fcn

    subroutine do_step()
       ii = ii + 1
       do jj = 1, nsteps
          do i = 2, s
             g(:,i) = g(:,1)
          enddo
          do j = 2, recur_steps
             g = matmul(g,a)
             do i = 1,s
                call fcn(y + g(:,i)*step, g(:,i))
             end do
          end do
          y = y + matmul(g,b)*step
          call fcn(y, g(:,1))
       enddo
       jlb_table(l)%data(0,ii) = y(2)
       jlb_table(l)%data(1,ii) = y(3)*scale1 
       jlb_table(l)%data(2,ii) = g(3,1)*scale2
     end subroutine do_step

  end subroutine jlb_tabulate

  subroutine jlb_get_jl(l, x, jl)
    integer l, i
    real(dl) x, jl, t
    if(l.lt. jlb_lmin)then
       call jlb_exact_lowl(l, x, jl)
    else
       t = (x-jlb_table(l)%xmin)/jlb_table(l)%dx
       i = floor(t)
       if(i.lt.0)then
          call jlb_small(l, x, jl) 
       elseif(i .ge. jlb_table(l)%n)then
          call jlb_large(l, x, jl)
       else
          t = t-i
          jl = -jlb_table(l)%data(0,i)*(-1+t)**3*(1+t*(3+6*t)) &
               +0.5d0*t*( &
               -2*jlb_table(l)%data(1,i)*(-1+t)**3*(1+3*t) &
               +t*(-jlb_table(l)%data(2,i)*(-1+t)**3+ &
               t*((jlb_table(l)%data(1,i+1)*(8-6*t) &
               +jlb_table(l)%data(2,i+1)*(-1+t))*(-1+t) &
               +2*jlb_table(l)%data(0,i+1)*(10+3*t*(-5+2*t)))) &
               )
       endif
    endif
  end subroutine jlb_get_jl

  subroutine jlb_get_jl_and_jlp(l, x, jl, jlp)
    integer l, i
    real(dl) x, jl, jlp, t
    if(l.lt. jlb_lmin)then
       call jlb_exact_lowl(l, x, jl, jlp)
    else
       t = (x-jlb_table(l)%xmin)/jlb_table(l)%dx
       i = floor(t)
       if(i.lt.0)then
          call jlb_small(l, x, jl, jlp) 
       elseif(i .ge. jlb_table(l)%n)then
          call jlb_large(l, x, jl, jlp)
       else
          t = t-i
          jl = -jlb_table(l)%data(0,i)*(-1+t)**3*(1+t*(3+6*t)) &
               +0.5d0*t*( &
               -2*jlb_table(l)%data(1,i)*(-1+t)**3*(1+3*t) &
               +t*(-jlb_table(l)%data(2,i)*(-1+t)**3+ &
               t*((jlb_table(l)%data(1,i+1)*(8-6*t) &
               +jlb_table(l)%data(2,i+1)*(-1+t))*(-1+t) &
               +2*jlb_table(l)%data(0,i+1)*(10+3*t*(-5+2*t)))) &
               )
          jlp = -jlb_table(l)%data(1,i)*(-1+t)**2*(-1+t*(-2+15*t)) &
               -0.5d0*t*(jlb_table(l)%data(2,i)*(-1+t)**2*(-2+5*t) &
               + t*(24*jlb_table(l)%data(1,i+1)-3*jlb_table(l)%data(2,i+1)+60*(jlb_table(l)%data(0,i)-jlb_table(l)%data(0,i+1))*(-1+t)**2 &
               + t*(-56*jlb_table(l)%data(1,i+1)+8*jlb_table(l)%data(2,i+1)+ t*(30*jlb_table(l)%data(1,i+1)-5*jlb_table(l)%data(2,i+1)))))
       endif
    endif
  end subroutine jlb_get_jl_and_jlp


  subroutine jlb_get_amp_phase(l, x, amp, phase)
    real(dl),parameter::  params(10) = (/ 0.230490766, -0.059341062, -0.008725434, 0.008935870, 0.277038862, -0.721654426, 1.311876363, 0.585017752, 2.675037879, -0.162880886 /)
    integer l
    real(dl) amp, phase
    real(dl) x,  nu2, nu, chb, ch2b, sh2b, shb, sx, sat, rat, z1, z2, fac, sqrtnu
    z1 = sphericalbesselj_zero(l, 1)
    if(x .le. z1)then
       amp = 0
       phase = coop_pio2
       return
    endif
    z2 = sphericalbesselj_zero(l, 2)
    nu2 = l*(l+1.d0)
    nu = sqrt(nu2)
    chb = x/nu
    ch2b = x**2/nu2
    sh2b = ch2b - 1.d0
    shb = sqrt(sh2b)
    sx = shb * nu
    sat = 1.d0/shb/sh2b/nu
    rat = sat*sat
    amp = dexp( &
         (1.d0/16-3.d0/8.d0*ch2b  &
         + ( 1.d0/128.d0 + ch2b*(3.d0/32.d0 + ch2b * (39/16.d0 + ch2b * (15.d0/8.d0))) &
         - (2523.d0/32.d0 + ch2b*(1521.d0/16.d0 + ch2b * (315.d0/16.d0)))*ch2b**3 * rat)*rat)*rat &
         -0.5d0*dlog(x*sx) &
         )
    if(x.lt. z2)then
       sqrtnu = sqrt(nu)
       fac = ( 1.d0 + tanh((z2-z1)*(exp(params(1))/(z1-x) + exp(params(2))/(z2-x)) + params(3) + params(5)/sqrtnu + params(7)/nu + params(10)/nu2) )/2.d0
       amp = amp * fac**exp(params(4)+params(6)/sqrtnu + params(8)/nu + params(9)/nu2)
    endif
    phase = (-coop_pio4 &
         + sx &
         -((-1.d0/128/nu2 + 0.125)/nu +nu)*acos(nu/x) &
         + sat*( &
         -1.d0/12.d0-ch2b/8.d0  &
         + rat * (-31/5760.d0 + ch2b *( 1062/5760.d0 + ch2b * (4644/5760.d0+ ch2b * ( -195/5760.d0 + 45/5760.d0 * ch2b))) &
         - rat * ch2b**2* (174126/32256.d0 + ch2b*(511445/32256.d0+ ch2b*(180144/32256.d0)) &
         ))))
  end subroutine jlb_get_amp_phase



end module coop_sbj_mod


