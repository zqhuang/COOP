module coop_jl_mod
  use coop_wrapper_typedef
  implicit none

#include "constants.h"

  private


  COOP_INT,parameter::coop_jl_lmax = 5101
  COOP_INT,parameter::coop_jl_lmin = 10

  !!coop_jl is much much more accurate than coop_SphericalBesselJ
  !!after first call of coop_jl (which requires initialization), coop_jl is also about 3 times faster than coop_sphericalbesselJ
  !! jl = coop_jl(l, x) ,   x can be real/double, or an array of real/double
  public::coop_jl, coop_jl_zero, coop_jljl_slow, coop_jl_startpoint, coop_jl_get_jl_and_jlp, coop_jl_check_init, coop_jl_get_jl, coop_jl_get_amp_phase, coop_jl_setup_amp_phase

  type coop_jl_table
     logical::tabulated = .false.
     COOP_REAL::xmax, xmin, dx
     COOP_INT::n
     COOP_INT::l
     COOP_INT::count = 0
     COOP_REAL,dimension(:,:),allocatable::data
   contains
     procedure::init => coop_jl_table_initialize
     procedure::free => coop_jl_table_free
  end type coop_jl_table
  
  type(coop_jl_table) :: coop_jl_global_table(coop_jl_lmin:coop_jl_lmax)

  type(coop_function),dimension(coop_jl_lmax)::coop_jl_amp, coop_jl_phase

!!yes -> 1.3 GB mem
!!no ->  0.4 GB mem
#define  COOP_JL_SUPER_PRECISION COOP_NO 

  interface coop_jl
     module procedure coop_jl_s, coop_jl_v, coop_jl_ss, coop_jl_vs
  end interface coop_jl


contains

  subroutine coop_jl_get_amp_phase(l, x, amp, phase)  !!assuming x>0
    COOP_INT l
    COOP_REAL x, amp, phase, lnxbynu
    if(l.eq.0)then
       amp = 1.d0/x
       phase = x - coop_pio2
    else
       lnxbynu = log(x/(l+0.5d0))
       amp = coop_jl_amp(l)%eval(lnxbynu)/x
       phase = coop_jl_phase(l)%eval(lnxbynu)
    endif
  end subroutine coop_jl_get_amp_phase

  subroutine coop_jl_setup_amp_phase(l)
    COOP_INT::l, n, i
    COOP_REAL,dimension(:),allocatable::x, amp, phase
    COOP_REAL::jl, dx, jltry, xmin, z2, xp1, xp2, z1, peak1, peak2, ampz2, nu, lnxbynumin, lnxbynumax
    if(l.le.0)return
    if(coop_jl_amp(l)%n .gt. 0)return  !!already initialized
    call coop_jl_check_init(l)
    xp1 = (l + 0.5d0)*(1.d0+0.80992d0/(l+0.5d0)**(2.d0/3.d0))
    jl = xp1*coop_jl(l, xp1)
    dx = 0.1d0/(l+0.5d0)
    do while(dx .gt. 1.d-7)
       jltry = (xp1+ dx)*coop_jl(l, xp1 + dx)
       do while(jltry .gt. jl)
          xp1 = xp1 + dx
          jl = jltry
          jltry = (xp1+ dx)*coop_jl(l, xp1 + dx)
       enddo
       dx = dx/2.d0
    enddo
    peak1 = jl
    nu = l + 0.5d0
    xp2 = (coop_jl_zero(l, 1) + coop_jl_zero(l, 2))/2.d0
    dx = 0.1d0/nu
    do while(dx .gt. 1.d-7)
       jltry = (xp2+ dx)*coop_jl(l, xp2 + dx)
       if(jltry .lt. jl)then
          do while(jltry .lt. jl)
             xp2 = xp2 + dx
             jl = jltry
             jltry = (xp2+ dx)*coop_jl(l, xp2 + dx)
          enddo
       else
          jltry = (xp2 - dx)*coop_jl(l, xp2-dx)
          do while(jltry .lt. jl)
             xp2 = xp2 - dx
             jl = jltry
             jltry = (xp2- dx)*coop_jl(l, xp2 - dx)
          enddo          
       endif
       dx = dx/2.d0
    enddo
    peak2 = -jl

    z2 = coop_jl_zero(l, 2)
    call coop_jl_get_amp_phase_asymptotic(l, z2, ampz2, jltry)
    ampz2 = ampz2*z2
    z1 = coop_jl_zero(l, 1)

    call coop_jl_startpoint(l, xmin)
    lnxbynumin = log(xmin/nu)
    lnxbynumax = log(max(3.d0, 5.d2/nu))

    n = ceiling((lnxbynumax-lnxbynumin)*max(l,100))

    allocate(x(n), amp(n), phase(n))

    call coop_set_uniform(n, x, lnxbynumin, lnxbynumax)
    x = exp(x)*nu

    do i = 1, n
       if(x(i) .le. z2)then
          if(x(i) .le. xp1)then
             amp(i) = peak1
             jl = coop_jl(l, x(i))*x(i)
             phase(i) = - acos(min(jl/peak1, 1.d0))
          elseif(x(i) .ge. xp2)then
             amp(i) = (peak2*(z2 - x(i)) + ampz2*(x(i)-xp2))/(z2-xp2)
             jl = coop_jl(l, x(i))*x(i)
             phase(i)  = coop_2pi - acos(max(jl/amp(i), -1.d0))
          else
             amp(i) = (peak1*(xp2-x(i))+peak2*(x(i) - xp1))/(xp2-xp1)
             jl = coop_jl(l, x(i))*x(i)
             phase(i) = acos(max(min(jl/amp(i), 1.d0), -1.d0))
          endif
       else
          call coop_jl_get_amp_phase_asymptotic(l, x(i), amp(i), phase(i))
          amp(i) = amp(i)*x(i)
       endif
    enddo
    call coop_jl_amp(l)%init(n = n, xmin = lnxbynumin, xmax = lnxbynumax, f = amp, check_boundary = .false., fleft = peak1, fright = 1.d0, slopeleft = 0.d0, sloperight = 0.d0)
    call coop_jl_phase(l)%init(n = n, xmin = lnxbynumin, xmax =  lnxbynumax, f = phase, check_boundary = .false., fleft = -coop_pio2, fright = phase(n), slopeleft = 0.d0, sloperight = 1.d0)
    deallocate(x, amp, phase)
    if(l.ge.coop_jl_lmin)call coop_jl_global_table(l)%free()
  end subroutine coop_jl_setup_amp_phase


  function coop_jl_phase_analytic(l, x) result(phase)
    COOP_INT l
    COOP_REAL phase
    COOP_REAL x,  nu2, nu, chb, ch2b, sh2b, shb, sx, sat, rat
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
  end function coop_jl_phase_analytic

  function coop_jl_zero(l, k) result(z)
    COOP_INT k, l
    COOP_REAL z, t, zmin, zmax, phasemax, phase, phasemid
    COOP_REAL,parameter::c1(12) =   (/ -0.245127043d+01, 0.112632772d+02, -0.207925322d+02, 0.198543806d+02, -0.110527677d+02, 0.406769876d+01, -0.109067653d+01, 0.206080657d+00, -0.260196064d-1 , 0.212880935d-2 , -0.101897908d-3 , 0.216576354d-5  /)
    COOP_REAL, parameter::c2(12) =  (/ 0.579779690d+01, -0.280853646d+02, 0.625315957d+02, -0.805135267d+02, 0.629961403d+02, -0.307235514d+02, 0.973806423d+01, -0.206907154d+01, 0.293355044d+00, -0.265358216d-1 , 0.138530622d-2 , -0.317522214d-4  /)
    if(l.eq.0)then
       z = coop_pi*k
       return
    endif
    select case(k)
    case(1)
       t = (l+0.5d0)**0.25d0
       z = l + 0.5d0 + 2.d0*(l+2.5d0)**0.323d0  + coop_polyvalue(12,c1, t)/(l+0.5d0)
    case(2)
       t = (l+0.5d0)**0.25d0
       z = l + 0.5d0 + 3.51d0*(l+5.d0)**0.323d0 + coop_polyvalue(12,c2, t)/(l+0.5d0)
    case default
       t = (l+0.5d0)**0.25d0
       zmin = l + 0.5d0 + 3.51d0*(l+5.d0)**0.323d0 + coop_polyvalue(12,c2, t)/(l+0.5d0)
       phase = coop_pi * (k-0.5d0) 
       zmax = max(coop_jl_zero_asym(l, k), zmin+coop_pi)
       phasemax = coop_jl_phase_analytic(l, zmax)
       do while(phasemax .lt. phase)
          zmax = zmax + max(phase - phasemax, coop_pi)
          phasemax = coop_jl_phase_analytic(l, zmax)
       enddo
       do while(zmax - zmin .gt. 1.d-6)
          z = (zmax+zmin)/2.d0
          phasemid = coop_jl_phase_analytic(l, z)
          if(phasemid .ge. phase)then
             zmax = z
          else
             zmin = z
          endif
       enddo
       z = (zmin+zmax)/2.d0
    end select
  end function coop_jl_zero


!!the k-th zero of spherical Bessel function j_l
  function coop_jl_zero_asym(l, k) result(z)
    COOP_REAL z, nusq, zsq
    COOP_INT l, k
    z = coop_pi * (l/2.d0+k)
    nusq = (l +0.5d0)**2
    zsq = z**2
    z = z - l*(l+1.d0)/(2.d0*z)*(1.d0+ ((-31.d0+28.d0*nusq) + &
         (3779.d0+ nusq*(-3928.d0 + 1328.d0*nusq)+(-6277237.d0+ nusq*(6342972.d0 + nusq*(-2461680.d0+nusq*444736.d0)))/(224.d0*zsq))/(40.d0*zsq))/(48.d0*zsq))
  end function coop_jl_zero_asym



!!start where j_l is ~ 10^{-8}
  subroutine coop_jl_startpoint(l, xstart, eps)
    COOP_REAL xstart, lnnu, nu
    COOP_REAL,optional::eps
    COOP_INT l
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
       xstart = xstart*(eps/1.d-8)**(1.d0/l)
    endif
  end subroutine coop_jl_startpoint



!!very accurate for x < l - (a few) * l^{1/3}
 subroutine coop_jl_small(l, x, jl, jlp)
    COOP_INT l
    COOP_REAL,intent(IN)::x
    COOP_REAL nu, nu2, x2, cosb, cos2b, sinb, sin2b, sin3b,  nusin3b, invnu3b, expterm
    COOP_REAL,intent(OUT):: jl
    COOP_REAL,intent(OUT),optional:: jlp
    if(l.lt.coop_jl_lmin)then
       if(present(jlp))then
          call coop_jl_exact_lowl(l, x, jl, jlp)
       else
          call coop_jl_exact_lowl(l, x, jl)
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
       expterm = dexp(-log_gamma(l+1.5)+l*(dlog(x)-coop_ln2)-coop_ln2)*coop_sqrtpi 
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
  end subroutine coop_jl_small



  subroutine coop_jl_large(l, x, jl, jlp)
    COOP_INT l
    COOP_REAL,intent(IN)::x
    COOP_REAL,intent(OUT)::jl
    COOP_REAL,intent(OUT),optional::jlp
    COOP_REAL sx, expterm, alpha, x2, nu, nu2, sh2b, shb, ch2b, chb, sat, rat, expderv, alphaderv, cschb, csch2b, invnu2, invnu
    if(l.lt.coop_jl_lmin)then
       call coop_jl_exact_lowl(l, x, jl, jlp)
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
  end subroutine coop_jl_large

  subroutine coop_jl_exact_lowl(l, x, jl, jlp)
    COOP_INT l
    COOP_REAL x, jl, v2, x2
    COOP_REAL,optional::jlp
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
          stop "coop_jl_exact_lowl called for high l"
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
          stop "coop_jl_exact_lowl called for high l"
    end select
  end subroutine coop_jl_exact_lowl

  subroutine coop_jl_check_init(l, save_memory)
    COOP_INT l, ell
    logical,optional::save_memory
    if(l.lt.coop_jl_lmin) return
    if(present(save_memory))then
       if(save_memory)then
          if(.not. coop_jl_global_table(l)%tabulated)then
             do ell=coop_jl_lmin, coop_jl_lmax
                if(coop_jl_global_table(ell)%tabulated)then
                   if(coop_jl_global_table(ell)%count .lt. -10000)then
                      call coop_jl_global_table(ell)%free()
                   else
                      coop_jl_global_table(ell)%count = coop_jl_global_table(ell)%count - 50
                   endif
                endif
             enddo
             call coop_jl_global_table(l)%init(l=l)
          endif
          return
       endif
    endif
    if(.not. coop_jl_global_table(l)%tabulated)call coop_jl_global_table(l)%init(l=l)
  end subroutine coop_jl_check_init

  subroutine coop_jl_get_jl(l, x, jl)
    COOP_INT, intent(IN)::l
    COOP_REAL,intent(IN)::x
    COOP_REAL,intent(OUT)::jl
    COOP_INT i, ell
    COOP_REAL t
    if(l.lt. coop_jl_lmin)then
       call coop_jl_exact_lowl(l, x, jl)
    else
       call coop_jl_check_init(l)
       t = (x-coop_jl_global_table(l)%xmin)/coop_jl_global_table(l)%dx
       i = floor(t)
       if(i.lt.0)then
          call coop_jl_small(l, x, jl) 
       elseif(i .ge. coop_jl_global_table(l)%n)then
          call coop_jl_large(l, x, jl)
       else
          t = t-i
          jl = -coop_jl_global_table(l)%data(0,i)*(-1+t)**3*(1+t*(3+6*t)) &
               +0.5d0*t*( &
               -2*coop_jl_global_table(l)%data(1,i)*(-1+t)**3*(1+3*t) &
               +t*(-coop_jl_global_table(l)%data(2,i)*(-1+t)**3+ &
               t*((coop_jl_global_table(l)%data(1,i+1)*(8-6*t) &
               +coop_jl_global_table(l)%data(2,i+1)*(-1+t))*(-1+t) &
               +2*coop_jl_global_table(l)%data(0,i+1)*(10+3*t*(-5+2*t)))) &
               )
       endif
    endif
  end subroutine coop_jl_get_jl

  subroutine coop_jl_get_jl_and_jlp(l, x, jl, jlp)
    COOP_INT l, i
    COOP_REAL x, jl, jlp, t
    if(l.lt. coop_jl_lmin)then
       call coop_jl_exact_lowl(l, x, jl, jlp)
    else
       call coop_jl_check_init(l)
       t = (x-coop_jl_global_table(l)%xmin)/coop_jl_global_table(l)%dx
       i = floor(t)
       if(i.lt.0)then
          call coop_jl_small(l, x, jl, jlp) 
       elseif(i .ge. coop_jl_global_table(l)%n)then
          call coop_jl_large(l, x, jl, jlp)
       else
          t = t-i
          jl = -coop_jl_global_table(l)%data(0,i)*(-1+t)**3*(1+t*(3+6*t)) &
               +0.5d0*t*( &
               -2*coop_jl_global_table(l)%data(1,i)*(-1+t)**3*(1+3*t) &
               +t*(-coop_jl_global_table(l)%data(2,i)*(-1+t)**3+ &
               t*((coop_jl_global_table(l)%data(1,i+1)*(8-6*t) &
               +coop_jl_global_table(l)%data(2,i+1)*(-1+t))*(-1+t) &
               +2*coop_jl_global_table(l)%data(0,i+1)*(10+3*t*(-5+2*t)))) &
               )
          jlp = -coop_jl_global_table(l)%data(1,i)*(-1+t)**2*(-1+t*(-2+15*t)) &
               -0.5d0*t*(coop_jl_global_table(l)%data(2,i)*(-1+t)**2*(-2+5*t) &
               + t*(24*coop_jl_global_table(l)%data(1,i+1)-3*coop_jl_global_table(l)%data(2,i+1)+60*(coop_jl_global_table(l)%data(0,i)-coop_jl_global_table(l)%data(0,i+1))*(-1+t)**2 &
               + t*(-56*coop_jl_global_table(l)%data(1,i+1)+8*coop_jl_global_table(l)%data(2,i+1)+ t*(30*coop_jl_global_table(l)%data(1,i+1)-5*coop_jl_global_table(l)%data(2,i+1)))))
       endif
    endif
  end subroutine coop_jl_get_jl_and_jlp


  subroutine coop_jl_get_amp_phase_asymptotic(l, x, amp, phase) !!assuming x>l+0.5
    COOP_INT l
    COOP_REAL amp, phase
    COOP_REAL x,  nu2, nu, chb, ch2b, sh2b, shb, sx, sat, rat, fac, sqrtnu
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
    phase = (-coop_pio4 &
         + sx &
         -((-1.d0/128/nu2 + 0.125)/nu +nu)*acos(nu/x) &
         + sat*( &
         -1.d0/12.d0-ch2b/8.d0  &
         + rat * (-31/5760.d0 + ch2b *( 1062/5760.d0 + ch2b * (4644/5760.d0+ ch2b * ( -195/5760.d0 + 45/5760.d0 * ch2b))) &
         - rat * ch2b**2* (174126/32256.d0 + ch2b*(511445/32256.d0+ ch2b*(180144/32256.d0)) &
         ))))
  end subroutine coop_jl_get_amp_phase_asymptotic

  function coop_jl_s(l, x) result(jl)
    COOP_REAL jl
    COOP_INT l
    COOP_REAL x
    call coop_jl_get_jl(l, x, jl)
  end function coop_jl_s

  function coop_jl_v(l, x) result(jl)
    COOP_INT l, i
    COOP_REAL x(:)
    COOP_REAL jl(size(x))
    call coop_jl_get_jl(l, x(1), jl(1))
    !$omp parallel do
    do i=2, size(x)
       call coop_jl_get_jl(l, x(i), jl(i))
    enddo
    !$omp end parallel do
  end function coop_jl_v



  function coop_jl_ss(l, x) result(jl)
    COOP_SINGLE jl
    COOP_INT l
    COOP_SINGLE x
    jl = coop_jl_s(l, dble(x))
  end function coop_jl_ss

  function coop_jl_vs(l, x) result(jl)
    COOP_INT l, i
    COOP_SINGLE x(:)
    COOP_SINGLE jl(size(x))
    jl =coop_jl_v(l, dble(x))
  end function coop_jl_vs



  subroutine coop_jl_table_initialize(this, l, max_x) 
    class(coop_jl_table)::this
    COOP_REAL,optional::max_x
    COOP_INT l, nsteps, ii , jj, i, j
    COOP_REAL step, y(3), lsq, scale1, scale2, rl
    COOP_INT, parameter :: s = 3
    COOP_REAL::step_max 
    COOP_INT recur_steps 
    COOP_REAL, parameter ::   b(s) = (/ 5.d0/18.d0, 4.d0/9.d0, 5.d0/18.d0/)
    ! Butcher tableau for 6th order Gauss-Legendre method
    COOP_REAL, parameter :: a(s,s) =  reshape( (/ &
         5.d0/36.d0, 2.d0/9.d0 - 1.d0/sqrt(15.d0), 5.d0/36.d0 - 0.5d0/sqrt(15.d0), &
         5.d0/36.d0 + sqrt(15.d0)/24.d0, 2.d0/9.d0, 5.d0/36.d0 - sqrt(15.d0)/24.d0, &
         5.d0/36.d0 + 0.5d0/sqrt(15.d0), 2.d0/9.d0 + 1.d0/sqrt(15.d0), 5.d0/36.d0 /), (/ s, s /) )
    COOP_REAL g(3,s)
    if(present(max_x))then
       if(this%tabulated .and. this%xmax .ge. max_x) return
    else
       if(this%tabulated) return
#if COOP_JL_SUPER_PRECISION
       this%xmax = min(max(l*2.8d0, 600.d0), l+2500.d0)
#else
       this%xmax = min(max(l*2.2d0, 300.d0), l+1000.d0)   
#endif
    endif
    if(l.lt.coop_jl_lmin)return
    call this%free()
    this%tabulated = .true.
#if COOP_JL_SUPER_PRECISION
    step_max = 0.02d0 + l/12900.d0
#else
    step_max = 0.031d0 + (l/6800.d0)
#endif
    recur_steps = 7 + floor(step_max/0.08)
    call coop_jl_startpoint(l, this%xmin, 1.d-24)
    this%n = max(min(16384, ceiling((this%xmax-this%xmin)/(step_max))),64)
    this%dx = (this%xmax-this%xmin)/this%n
    scale1 = this%dx
    scale2 = scale1**2
    allocate(this%data(0:2, 0:this%n))
    lsq = l*(l+1.d0)
    y(1) = this%xmin
    call coop_jl_small(l, y(1), y(2), y(3))
    call fcn(y, g(:,1))
    ii = 0
    this%data(0,ii) = y(2)
    this%data(1,ii) = y(3)*scale1 
    this%data(2,ii) = g(3,1)*scale2

    nsteps = ceiling(this%dx/ min((this%xmin+ii*this%dx)/l, 1.d0)/min(step_max, 0.15d0))
    step = this%dx/nsteps
    rl = l*1.3d0
    do while(ii .lt. this%n)
       call do_step()
       if(y(1) .ge. rl )exit
       nsteps = ceiling(this%dx/(y(1)/rl)/min(step_max,0.15d0))
       step = this%dx/nsteps
    enddo

    do while(ii .lt. this%n)
       call do_step()
    enddo

  contains

    subroutine fcn(yy, yyp)
      COOP_REAL yy(3), yyp(3)
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
       this%data(0,ii) = y(2)
       this%data(1,ii) = y(3)*scale1 
       this%data(2,ii) = g(3,1)*scale2
     end subroutine do_step

  end subroutine coop_jl_table_initialize

  subroutine coop_jl_table_free(this)
    class(coop_jl_table)::this
    if(allocated(this%data))deallocate(this%data)
    this%tabulated = .false.
    this%count = 0
  end subroutine coop_jl_table_free



  function coop_jljl_slow(l, x1, x2) result(jljl)
    COOP_INT l
    COOP_REAL x1, x2, jljl
    COOP_REAL ph1, ph2, a1, a2
    call coop_jl_get_amp_phase_asymptotic(l, x1, a1, ph1)
    call coop_jl_get_amp_phase_asymptotic(l, x2, a2, ph2)
    jljl = a1*a2*cos(ph1-ph2)/2.d0
  end function coop_jljl_slow


end module coop_jl_mod


