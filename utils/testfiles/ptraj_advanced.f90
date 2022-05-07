!!in linux environment compile the code with
!! gfortran -O3 particle_trajectory.f90 -o calcp
!!and run it with
!! ./calcp > output.txt
!! output: t, x, y, z, vx, vy vz (units: microsecond, mm, mm, mm, m/s, m/s, m/s)
module pt_subs
  implicit none
  integer,parameter::IB = selected_int_kind(8)
  integer,parameter::DL = selected_real_kind(12)
  real(kind=DL),parameter::pi  = 3.14159265358979323846264338327950288_dl
  !!---------------------USER DEFINED VARIABLES ------------------------
  real(kind=DL),parameter::mass = 171.0_dl * 1.67e-27_dl !!SI unit
  real(kind=DL),parameter::charge = 1.6e-19_dl !!SI unit
  real(kind=DL),parameter::norm = 43503.0_dl 
  integer(kind=IB),parameter::nx(3) = (/ 101, 101, 126 /)
  real(kind=DL),parameter::xmin(3) = (/ -1.e-3_dl, -1.e-3_dl, 1.e-3_dl /) !!SI unit
  real(kind=DL),parameter::xmax(3) = (/ 1.e-3_dl, 1.e-3_dl, 3.5e-3_dl /) !!SI unit
  real(kind=DL),dimension(3),parameter::xini = (/ 0.0_dl, 0.0_dl, 2.1e-3_dl /) !!SI unit
  real(kind=DL),dimension(3),parameter::vini = (/ 1.0_dl, 1.0_dl, 1.0_dl /) !!SI unit
  real(kind=DL),parameter::freq = 2.e7_dl
  real(kind=DL),parameter::tend = 10._dl * 1.e-6_dl 
  !!--------------------------------------------------------------------  
  real(kind=DL),parameter::charge_by_mass = charge/mass  
  real(kind=DL),parameter::omega = 2.0_dl*pi*freq
  real(kind=DL)::E_field(3, nx(1), nx(2), nx(3)), potential(nx(1), nx(2), nx(3))
  real(kind=DL),parameter::dx(3) = (xmax-xmin)/(nx-1)
  
contains
  
  subroutine load_data(filename)
    character(LEN=*)::filename
    integer(kind=IB)::i,j,k
    character(LEN=1024)::str
    real(kind=DL) line(7)
    open(10, FILE=filename)
    read(10,"(a)") str
    print*, trim(str)    
    do while(str(1:5) .ne. "x,y,z")
       read(10,"(A)") str
       !       print*, trim(str)
    enddo
    do i=1, nx(1)
       do j=1, nx(2)
          do k=1, nx(3)
             read(10,*, END=100, ERR=100) line
             line(1:3) = line(1:3)*1.e-3_dl
             if(nint((line(1)-xmin(1))/dx(1))+1 .ne. i &
                  .or. nint((line(2)-xmin(2))/dx(2))+1 .ne. j &
                  .or. nint((line(3)-xmin(3))/dx(3))+1 .ne. k )then
                goto 100
             endif
             E_field(:, i, j, k) = line(4:6)
             potential(i, j, k) = line(7)
          enddo
       enddo
    enddo
    close(10)
    E_field = E_field * (norm * charge_by_mass*1.e3_dl) 
    potential = potential * (norm * charge_by_mass)
    open(10, FILE="saved.dat", FORM="unformatted")
    write(10) E_field
    write(10)potential
    close(10)
    return
100 stop "Error in the data file"
  end subroutine load_data


    subroutine load_data_quick()
    open(10, FILE="saved.dat",FORM="unformatted")
    read(10) E_field
    read(10) potential
    close(10)
  end subroutine load_data_quick


  subroutine get_acceleration_with_E(x, a)
    real(kind=DL)::x(3), a(3), u(3), d(3)
    integer(kind=IB)::i(3)
    u = (x - xmin)/dx+1.0_dl
    i = floor(u)
    if(any(i.le.0) .or. any(i.ge. nx))then
       write(*, *) x
       stop "Error: particle hits the domain boundary."
    endif
    u = u - i
    d = 1.0_dl - u
    !!use E, biliner interpolation    
    a = d(1)*( &
         d(2)*( &
             E_field(:, i(1)  , i(2)  , i(3)  ) * d(3) &
           + E_field(:, i(1)  , i(2)  , i(3)+1) * u(3) &
         ) &
         + u(2)*(  &
         E_field(:, i(1)  , i(2)+1, i(3)  ) * d(3) &
         + E_field(:, i(1)  , i(2)+1, i(3)+1) * u(3) &
         ) &
         ) &
         + u(1) * ( &
         d(2)*( &
         E_field(:, i(1)+1, i(2)  , i(3)  ) * d(3) &
         + E_field(:, i(1)+1, i(2)  , i(3)+1) * u(3) &
         ) &
         + u(2)*(  &
         E_field(:, i(1)+1, i(2)+1, i(3)  ) * d(3) &
         + E_field(:, i(1)+1, i(2)+1, i(3)+1) * u(3) &
         ) &
         ) 
  end subroutine get_acceleration_with_E


  !!use cubic interpolation in z direction only (motion in x, y directions is negligible anyway)
  subroutine get_acceleration_with_E_interpz3(x, a)
    real(kind=DL)::x(3), a(3), u(3), d(3)
    integer(kind=IB)::i(3)
    u = (x - xmin)/dx+1.0_dl
    i = floor(u)
    if(any(i.le.1) .or. any(i.ge. nx-1))then
       write(*, *) x
       stop "Error: particle hits the domain boundary."
    endif
    u = u - i
    d = 1.0_dl - u
    a = (1.0_dl + u(3))*(1._dl - u(3)/2._dl)* d(3) *( &
         d(2)*( &
             E_field(:, i(1), i(2), i(3)) * d(1) &
           + E_field(:, i(1)+1, i(2), i(3)) * u(1) &
         ) &
         + u(2)*(  &
         E_field(:, i(1), i(2)+1, i(3)) * d(1) &
         + E_field(:, i(1)+1, i(2)+1, i(3)) * u(1) &
         ) &
         ) &
         + u(3)*(1.0_dl + u(3) * d(3) / 2._dl ) * ( &
         d(2)*( &
             E_field(:, i(1), i(2), i(3)+1) * d(1) &
           + E_field(:, i(1)+1, i(2), i(3)+1) * u(1) &
         ) &
         + u(2)*(  &
         E_field(:, i(1), i(2)+1, i(3)+1) * d(1) &
         + E_field(:, i(1)+1, i(2)+1, i(3)+1) * u(1) &
         ) &
         ) &
         + u(3)*(-1.0_dl/3.0_dl + u(3) *(0.5_dl -  u(3)/ 6._dl)) * ( &
         d(2)*( &
             E_field(:, i(1), i(2), i(3)-1) * d(1) &
           + E_field(:, i(1)+1, i(2), i(3)-1) * u(1) &
         ) &
         + u(2)*(  &
         E_field(:, i(1), i(2)+1, i(3)-1) * d(1) &
         + E_field(:, i(1)+1, i(2)+1, i(3)-1) * u(1) &
         ) &
         ) &
         + u(3)/6._dl * (u(3)**2 - 1.0_dl)  * ( &
         d(2)*( &
             E_field(:, i(1), i(2), i(3)+2) * d(1) &
           + E_field(:, i(1)+1, i(2), i(3)+2) * u(1) &
         ) &
         + u(2)*(  &
         E_field(:, i(1), i(2)+1, i(3)+2) * d(1) &
         + E_field(:, i(1)+1, i(2)+1, i(3)+2) * u(1) &
         ) &
         ) 
  end subroutine get_acceleration_with_E_interpz3
  

  !!use cubic interpolation in z direction only (motion in x, y directions is negligible anyway)
  subroutine get_acceleration_with_E_interp3(x, a)
    real(kind=DL)::x(3), a(3), r(3), w(4, 3)
    integer(kind=IB)::i(3), j
    r = (x - xmin)/dx+1.0_dl
    i = floor(r)
    if(any(i.le.1) .or. any(i.ge. nx-1))then
       write(*, *) x
       stop "Error: particle hits the domain boundary."
    endif
    r = r - i
    w(1, :) = r*(-1._dl/3._dl + r*(0.5_dl - r/6._dl))
    w(2, :) = (1._dl - r**2)*(1._dl - r/2._dl)
    w(3, :) = r*(1._dl + r*(1._dl - r)/2._dl)
    w(4, :) = r*(r**2-1._dl)/6._dl
    !$omp parallel do
    do j=1, 3
       a(j) = dot_product(w(:, 1), &
            matmul(&
            (&
            E_field(j, i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3)-1)*w(1, 3) &
            + E_field(j, i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3))*w(2, 3) &
            + E_field(j, i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3)+1)*w(3, 3) &
            + E_field(j, i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3)+2)*w(4, 3) &       
            )  &
            ,w(:, 2)) &
            )
    enddo
    !$omp end parallel do
  end subroutine get_acceleration_with_E_interp3


  !!interp3 algorithm
  !! f = f(-1) * (   -r/3 + r^2/2 - r^3/6) &
  !!   + f(0)  * (1 - r/2 - r^2   + r^3/2) &
  !!   + f(1)  * (    r   + r^2/2 - r^3/2) &
  !!   + f(2)  * (   -r/6         + r^3/6)
  !!f' = f(-1) * (  -1/3 + r  - r^2/2) &
  !!   + f(0)  * ( - 1/2 - 2r + 3r^2/2) &
  !!   + f(1)  * (     1 + r  - 3r^2/2) &
  !!   + f(2)  * (  -1/6      + r^2/2)
  function potential_interp(x)
    real(kind=DL)::x(3), potential_interp, r(3), w(4, 3)
    integer(kind=IB)::i(3)
    r = (x - xmin)/dx+1.0_dl
    i = floor(r)
    if(any(i.le.1) .or. any(i.ge. nx-1))then
       potential = 1.e99_dl
       return
    endif
    r =r - i
    w(1, :) = r*(-1._dl/3._dl + r*(0.5_dl - r/6._dl))
    w(2, :) = (1._dl - r**2)*(1._dl - r/2._dl)
    w(3, :) = r*(1._dl + r*(1._dl - r)/2._dl)
    w(4, :) = r*(r**2-1._dl)/6._dl
    potential_interp = dot_product(w(:, 1), &
         matmul(&
         (&
         potential(i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3)-1)*w(1, 3) &
         + potential(i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3))*w(2, 3) &
         + potential(i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3)+1)*w(3, 3) &
         + potential(i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3)+2)*w(4, 3) &       
         )  &
         ,w(:, 2)) &
         )
  end function potential_interp
  
  !!interp3 algorithm
  !! f = f(-1) * (   -r/3 + r^2/2 - r^3/6) &
  !!   + f(0)  * (1 - r/2 - r^2   + r^3/2) &
  !!   + f(1)  * (    r   + r^2/2 - r^3/2) &
  !!   + f(2)  * (   -r/6         + r^3/6)
  !!f' = f(-1) * (  -1/3 + r  - r^2/2) &
  !!   + f(0)  * ( - 1/2 - 2r + 3r^2/2) &
  !!   + f(1)  * (     1 + r  - 3r^2/2) &
  !!   + f(2)  * (  -1/6      + r^2/2)
  subroutine get_acceleration_with_U(x, a)
    real(kind=DL)::x(3), a(3), r(3), w(4, 3), wp(4, 3)
    integer(kind=IB)::i(3)
    r = (x - xmin)/dx+1.0_dl
    i = floor(r)
    if(any(i.le.1) .or. any(i.ge. nx-1))then
       write(*, *) x
       stop "Error: particle hits the domain boundary."
    endif
    r =r - i
    w(1, :) = r*(-1._dl/3._dl + r*(0.5_dl - r/6._dl))
    w(2, :) = (1._dl - r**2)*(1._dl - r/2._dl)
    w(3, :) = r*(1._dl + r*(1._dl - r)/2._dl)
    w(4, :) = r*(r**2-1._dl)/6._dl
    wp(1, :) = -(-1._dl/3._dl + r*(1._dl - r/2._dl))/dx
    wp(2, :) = -(-0.5_dl + r*(-2._dl + r * 1.5_dl))/dx
    wp(3, :) = -(1._dl + r * ( 1._dl - 1.5_dl * r))/dx
    wp(4, :) = -(-1._dl/6._dl + r**2/2._dl)/dx
    a(1) = dot_product(wp(:, 1), &
         matmul(&
         (&
         potential(i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3)-1)*w(1, 3) &
         + potential(i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3))*w(2, 3) &
         + potential(i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3)+1)*w(3, 3) &
         + potential(i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3)+2)*w(4, 3) &       
         )  &
         ,w(:, 2)) &
         )
    a(2) = dot_product(w(:, 1), &
         matmul(&
         (&
         potential(i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3)-1)*w(1, 3) &
         + potential(i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3))*w(2, 3) &
         + potential(i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3)+1)*w(3, 3) &
         + potential(i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3)+2)*w(4, 3) &       
         )  &
         ,wp(:, 2)) &
         )
    a(3) = dot_product(w(:, 1), &
         matmul(&
         (&
         potential(i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3)-1)*wp(1, 3) &
         + potential(i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3))*wp(2, 3) &
         + potential(i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3)+1)*wp(3, 3) &
         + potential(i(1)-1:i(1)+2, i(2)-1:i(2)+2, i(3)+2)*wp(4, 3) &       
         )  &
         ,w(:, 2)) &
         )
  end subroutine get_acceleration_with_U


  !!y(1:3) particle position
  !!y(4:6) particle velocity
  subroutine eom(n, t, y, yp)
    integer(kind=IB)::n
    real(kind=DL)::t, y(n), yp(n)
    yp(1:3) = y(4:6)
!!    call get_acceleration_with_E_interp3(y(1:3), yp(4:6))
    call get_acceleration_with_E(y(1:3), yp(4:6))
!!    call get_acceleration_with_U(y(1:3), yp(4:6))        
    yp(4:6) = yp(4:6) *sin(omega*t) 
  end subroutine eom


  !! given fcn(n, t, y(1:n), yp(1:n)) (input n,t, y, return yp = dy/dt);
  !! evolve t to t+h and update y
  subroutine RungeKutta4th(n, fcn, t, y, h)
    external fcn
    integer(kind=IB)  n
    real(kind=DL) y(n)
    real(kind=DL) t, h
    real(kind=DL) k1(n), k2(n), k3(n), k4(n)
    call fcn(n, t, y, k1)
    call fcn(n, t+h/2._dl, y+k1*(h/2._dl), k2)
    call fcn(n, t+h/2._dl, y+k2*(h/2._dl), k3)
    t = t + h
    call fcn(n, t, y+k3*h, k4)
    y = y+(k1+2._dl*(k2+k3)+k4)*(h/6._dl)
  end subroutine RungeKutta4th
  

  !! given fcn(n, t, y(1:n), yp(1:n)) (input n,t, y, return yp = dy/dt);
  !! evolve t to t+h and update y  
  subroutine RungeKutta6th(n, fcn, t, y, h)
    external fcn
    !! fcn(n, t, y(1:n), yp(1:n))
    real(kind=DL),parameter::sqrt21 = sqrt(21._dl)
    integer(kind=IB)  n
    real(kind=DL) y(n)
    real(kind=DL) t, h
    real(kind=DL) k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), k7(n)
    call fcn(n, t, y, k1)  !!0
    call fcn(n, t+h, y+k1*h, k2)  !! 1
    call fcn(n, t+h/2._dl, y+(3._dl*k1+k2)*(h/8._dl), k3) !!0.5
    call fcn(n, t+(2._dl/3._dl)*h, y+ (4._dl*(k1+k3)+k2)*((2._dl/27._dl)*h), k4) !! 0.67
    call fcn(n, t+((7._dl-sqrt21)/14._dl)*h, &
         y+((3._dl*(3._dl*sqrt21 - 7._dl)/392._dl)*k1 - (8._dl*(7._dl-sqrt21)/392._dl)*k2 &
         + (48._dl*(7._dl-sqrt21)/392._dl)*k3 - (3._dl*(21._dl - sqrt21)/392._dl)*k4)*h, k5)  !! 0.17
    call fcn(n, t+((7._dl+sqrt21)/14._dl)*h, &
         y+ ((-5._dl*(231._dl+51._dl*sqrt21)/1960._dl)*k1 &
         - (40._dl*(7._dl+sqrt21)/1960._dl)*k2 - (320._dl*sqrt21/1960._dl)*k3 &
         + (3._dl*(21._dl + 121._dl*sqrt21)/1960._dl)*k4 + (392._dl*(6._dl+sqrt21)/1960._dl)*k5)*h, k6)  !!0.8
    t = t+h
    call fcn(n, t, y+(((22._dl+7._dl*sqrt21)/12._dl)*k1 &
         +(2._dl/3._dl)*k2 + (2._dl*(7._dl*sqrt21-5._dl)/9._dl)*k3 &
         - (7._dl*(3._dl*sqrt21-2._dl)/20._dl)*k4 - (7._dl*(49._dl+9._dl*sqrt21)/90._dl)*k5 &
         + (7._dl*(7._dl-sqrt21)/18._dl)*k6)*h, k7) !! 1
    y = y+  ((k1+k7)/20._dl + (49._dl/180._dl)*(k5+k6) + (16._dl/45._dl)*k3)*h
  end subroutine RungeKutta6th


  !! given fcn(n, t, y(1:n), yp(1:n)) (input n,t, y, return yp = dy/dt);
  !! evolve t to t+h and update y  
  subroutine RungeKutta8th (n, fcn, t, y, h)
    real(kind=DL),parameter::sqrt6 = sqrt(6.0_dl)
    real(kind=DL),parameter:: c1 = 103._dl / 1680._dl
    real(kind=DL),parameter:: c8 = -27._dl / 140._dl
    real(kind=DL),parameter:: c9 = 76._dl / 105._dl
    real(kind=DL),parameter:: c10 = -201._dl / 280._dl
    real(kind=DL),parameter:: c11 = 1024._dl / 1365._dl
    real(kind=DL),parameter:: c12 = 3._dl / 7280._dl
    real(kind=DL),parameter:: c13 = 12._dl / 35._dl
    real(kind=DL),parameter:: c14 = 9._dl / 280._dl
    real(kind=DL),parameter:: a2 = 1._dl / 12._dl
    real(kind=DL),parameter:: a3 = 1._dl / 9._dl
    real(kind=DL),parameter:: a4 = 1._dl / 6._dl
    real(kind=DL),parameter:: a5 = 2._dl * (1._dl + sqrt6 ) / 15._dl
    real(kind=DL),parameter:: a6 = (6._dl + sqrt6 ) / 15._dl
    real(kind=DL),parameter:: a7 = (6._dl - sqrt6 ) / 15._dl
    real(kind=DL),parameter:: a8 = 2._dl / 3._dl
    real(kind=DL),parameter:: a9 = 1._dl / 2._dl
    real(kind=DL),parameter:: a10 = 1._dl / 3._dl
    real(kind=DL),parameter:: a11 = 1._dl / 4._dl
    real(kind=DL),parameter:: a12 = 4._dl / 3._dl
    real(kind=DL),parameter:: a13 = 5._dl / 6._dl
    real(kind=DL),parameter:: a15 = 1._dl / 6._dl
    real(kind=DL),parameter:: b31 = 1._dl / 27._dl
    real(kind=DL),parameter:: b32 = 2._dl / 27._dl
    real(kind=DL),parameter:: b41 = 1._dl / 24._dl
    real(kind=DL),parameter:: b43 = 3._dl / 24._dl
    real(kind=DL),parameter:: b51 = (4._dl + 94._dl * sqrt6 ) / 375._dl
    real(kind=DL),parameter:: b53 = -(282._dl + 252._dl * sqrt6 ) / 375._dl
    real(kind=DL),parameter:: b54 = (328._dl + 208._dl * sqrt6 ) / 375._dl
    real(kind=DL),parameter:: b61 = (9._dl - sqrt6 ) / 150._dl
    real(kind=DL),parameter:: b64 = (312._dl + 32._dl * sqrt6 ) / 1425._dl
    real(kind=DL),parameter:: b65 = (69._dl + 29._dl * sqrt6 ) / 570._dl
    real(kind=DL),parameter:: b71 = (927._dl - 347._dl * sqrt6 ) / 1250._dl
    real(kind=DL),parameter:: b74 = (-16248._dl + 7328._dl * sqrt6 ) / 9375._dl
    real(kind=DL),parameter:: b75 = (-489._dl + 179._dl * sqrt6 ) / 3750._dl
    real(kind=DL),parameter:: b76 = (14268._dl - 5798._dl * sqrt6 ) / 9375._dl
    real(kind=DL),parameter:: b81 = 4._dl / 54._dl
    real(kind=DL),parameter:: b86 = (16._dl - sqrt6 ) / 54._dl
    real(kind=DL),parameter:: b87 = (16._dl + sqrt6 ) / 54._dl
    real(kind=DL),parameter:: b91 = 38._dl / 512._dl
    real(kind=DL),parameter:: b96 = (118._dl - 23._dl * sqrt6 ) / 512._dl
    real(kind=DL),parameter:: b97 = (118._dl + 23._dl * sqrt6 ) / 512._dl
    real(kind=DL),parameter:: b98 = - 18._dl / 512._dl
    real(kind=DL),parameter:: b10_1 = 11._dl / 144._dl
    real(kind=DL),parameter:: b10_6 = (266._dl - sqrt6 ) / 864._dl
    real(kind=DL),parameter:: b10_7 = (266._dl + sqrt6 ) / 864._dl
    real(kind=DL),parameter:: b10_8 = - 1._dl / 16._dl
    real(kind=DL),parameter:: b10_9 = - 8._dl / 27._dl
    real(kind=DL),parameter:: b11_1 = (5034._dl - 271._dl * sqrt6 ) / 61440._dl
    real(kind=DL),parameter:: b11_7 = (7859._dl - 1626._dl * sqrt6 ) / 10240._dl
    real(kind=DL),parameter:: b11_8 = (-2232._dl + 813._dl * sqrt6 ) / 20480._dl
    real(kind=DL),parameter:: b11_9 = (-594._dl  + 271._dl * sqrt6 ) / 960._dl
    real(kind=DL),parameter:: b11_10 = (657._dl - 813._dl * sqrt6 ) / 5120._dl
    real(kind=DL),parameter:: b12_1 = (5996._dl - 3794._dl * sqrt6 ) / 405._dl
    real(kind=DL),parameter:: b12_6 = (-4342._dl - 338._dl * sqrt6 ) / 9._dl
    real(kind=DL),parameter:: b12_7 = (154922._dl - 40458._dl * sqrt6 ) / 135._dl
    real(kind=DL),parameter:: b12_8 = (-4176._dl + 3794._dl * sqrt6 ) / 45._dl
    real(kind=DL),parameter:: b12_9 = (-340864._dl + 242816._dl * sqrt6 ) / 405._dl
    real(kind=DL),parameter:: b12_10 = (26304._dl - 15176._dl * sqrt6 ) / 45._dl
    real(kind=DL),parameter:: b12_11 = -26624._dl / 81._dl
    real(kind=DL),parameter:: b13_1 = (3793._dl + 2168._dl * sqrt6 ) / 103680._dl
    real(kind=DL),parameter:: b13_6 = (4042._dl + 2263._dl * sqrt6 ) / 13824._dl
    real(kind=DL),parameter:: b13_7 = (-231278._dl + 40717._dl * sqrt6 ) / 69120._dl
    real(kind=DL),parameter:: b13_8 = (7947._dl - 2168._dl * sqrt6 ) / 11520._dl
    real(kind=DL),parameter:: b13_9 = (1048._dl - 542._dl * sqrt6 ) / 405._dl
    real(kind=DL),parameter:: b13_10 = (-1383._dl + 542._dl * sqrt6 ) / 720._dl
    real(kind=DL),parameter:: b13_11 = 2624._dl / 1053._dl
    real(kind=DL),parameter:: b13_12 = 3._dl / 1664._dl
    real(kind=DL),parameter:: b14_1 = -137._dl / 1296._dl
    real(kind=DL),parameter:: b14_6 = (5642._dl - 337._dl * sqrt6 ) / 864._dl
    real(kind=DL),parameter:: b14_7 = (5642._dl + 337._dl * sqrt6 ) / 864._dl
    real(kind=DL),parameter:: b14_8 = -299._dl / 48._dl
    real(kind=DL),parameter:: b14_9 = 184._dl / 81._dl
    real(kind=DL),parameter:: b14_10 = -44._dl / 9._dl
    real(kind=DL),parameter:: b14_11 = -5120._dl / 1053._dl
    real(kind=DL),parameter:: b14_12 = -11._dl / 468._dl
    real(kind=DL),parameter:: b14_13 = 16._dl / 9._dl
    real(kind=DL),parameter:: b15_1 = (33617._dl - 2168._dl * sqrt6 ) / 518400._dl
    real(kind=DL),parameter:: b15_6 = (-3846._dl + 31._dl * sqrt6 ) / 13824._dl
    real(kind=DL),parameter:: b15_7 = (155338._dl - 52807._dl * sqrt6 ) / 345600._dl
    real(kind=DL),parameter:: b15_8 = (-12537._dl + 2168._dl * sqrt6 ) / 57600._dl
    real(kind=DL),parameter:: b15_9 = (92._dl + 542._dl * sqrt6 ) / 2025._dl
    real(kind=DL),parameter:: b15_10 = (-1797._dl - 542._dl * sqrt6 ) / 3600._dl
    real(kind=DL),parameter:: b15_11 = 320._dl / 567._dl
    real(kind=DL),parameter:: b15_12 = -1._dl / 1920._dl
    real(kind=DL),parameter:: b15_13 = 4._dl / 105._dl
    real(kind=DL),parameter:: b16_1 = (-36487._dl - 30352._dl * sqrt6 ) / 279600._dl
    real(kind=DL),parameter:: b16_6 = (-29666._dl - 4499._dl * sqrt6 ) / 7456._dl
    real(kind=DL),parameter:: b16_7 = (2779182._dl - 615973._dl * sqrt6 ) / 186400._dl
    real(kind=DL),parameter:: b16_8 = (-94329._dl + 91056._dl * sqrt6 ) / 93200._dl
    real(kind=DL),parameter:: b16_9 = (-232192._dl + 121408._dl * sqrt6 ) / 17475._dl
    real(kind=DL),parameter:: b16_10 = (101226._dl - 22764._dl * sqrt6 ) / 5825._dl
    real(kind=DL),parameter:: b16_11 = - 169984._dl / 9087._dl
    real(kind=DL),parameter:: b16_12 = - 87._dl / 30290._dl
    real(kind=DL),parameter:: b16_13 =  492._dl / 1165._dl
    real(kind=DL),parameter:: b16_15 =  1260._dl / 233._dl
    real(kind=DL),parameter:: e1 = -1911._dl / 109200._dl
    real(kind=DL),parameter:: e8 = 34398._dl / 109200._dl
    real(kind=DL),parameter:: e9 = -61152._dl / 109200._dl
    real(kind=DL),parameter:: e10 = 114660._dl / 109200._dl
    real(kind=DL),parameter:: e11 = -114688._dl / 109200._dl
    real(kind=DL),parameter:: e12 = -63._dl / 109200._dl
    real(kind=DL),parameter:: e13 = -13104._dl / 109200._dl
    real(kind=DL),parameter:: e14 = -3510._dl / 109200._dl
    real(kind=DL),parameter:: e15 = 39312._dl / 109200._dl
    real(kind=DL),parameter:: e16 = 6058._dl / 109200._dl
    integer(kind=IB)  n
    external fcn
    !! fcn(n, t, y(1:n), yp(1:n))
    real(kind=DL) t, y(n), h
    real(kind=DL),dimension(n)::k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16
    real(kind=DL) h6, h12
    h12 = a2 * h
    h6 = a3 * h
    call fcn(n, t, y, k1)
    call fcn(n, t+h12, y +  h12 * k1, k2)
    call fcn(n, t+a3*h, y +  h * ( b31*k1 + b32*k2), k3 )
    call fcn(n, t+a4*h, y +  h * ( b41*k1 + b43*k3), k4 )
    call fcn(n, t+a5*h, y +  h * ( b51*k1 + b53*k3 + b54*k4), k5 )
    call fcn(n, t+a6*h, y +  h * ( b61*k1 + b64*k4 + b65*k5), k6 )
    call fcn(n, t+a7*h, y +  h * ( b71*k1 + b74*k4 + b75*k5 + b76*k6), k7)
    call fcn(n, t+a8*h, y +  h * ( b81*k1 + b86*k6 + b87*k7), k8 )
    call fcn(n, t+a9*h, y +  h * ( b91*k1 + b96*k6 + b97*k7 + b98*k8), k9 )
    call fcn(n, t+a10*h, y +  h * ( b10_1*k1 + b10_6*k6 + b10_7*k7 + b10_8*k8 &
         + b10_9*k9 ), k10 )
    call fcn(n, t+a11*h, y +  h * ( b11_1*k1 + b11_7*k7 + b11_8*k8 + b11_9*k9 &
         + b11_10 * k10 ), k11 )
    call fcn(n, t+a12*h, y +  h * ( b12_1*k1 + b12_6*k6 + b12_7*k7 + b12_8*k8 &
         + b12_9*k9 + b12_10 * k10 + b12_11 * k11 ), k12 )
    call fcn(n, t+a13*h, y +  h * ( b13_1*k1 + b13_6*k6 + b13_7*k7 + b13_8*k8 &
         + b13_9*k9 + b13_10*k10 + b13_11*k11 + b13_12*k12 ), k13 )
    call fcn(n, t+h, y +  h * ( b14_1*k1 + b14_6*k6 + b14_7*k7 + b14_8*k8 &
         + b14_9*k9 + b14_10*k10 + b14_11*k11 + b14_12*k12 + b14_13*k13 ), k14 )
   call fcn(n, t+a15*h, y +  h * ( b15_1*k1 + b15_6*k6 + b15_7*k7 + b15_8*k8 &
        + b15_9*k9 + b15_10*k10 + b15_11*k11 + b15_12*k12 + b15_13*k13 ), k15 )
   t = t + h
   call fcn(n, t, y +  h * ( b16_1*k1 + b16_6*k6 + b16_7*k7 + b16_8*k8 &
        + b16_9*k9 + b16_10*k10 + b16_11*k11 + b16_12*k12 + b16_13*k13 &
        + b16_15*k15), k16 )
   y = y +  h * ( c1 * k1 + c8 * k8 + c9 * k9 + c10 * k10 + c11 * k11 &
        + c12 * k12 + c13 * k13 + c14 * k14 )
 end subroutine RungeKutta8th

end module pt_subs


program pt
  use pt_subs
  integer(kind=IB),parameter::n=6
  real(kind=DL)::y(n), t, yp(n)
  integer(kind=IB)::i
  real(kind=DL),parameter::stepsize = 2.e-4_dl/freq
  !call load_data("19_05_18_01.csv")
  call load_data_quick()
  y(1:3) = xini
  y(4:6) = vini
  t = 0.0_dl
  write(*, "(7G14.5)") t*1.e6, y(1:3)*1.e3, y(4:6) !!output unit: microsecond, mm, m/s !!test
  do i=1, nint(tend/stepsize)
     call RungeKutta4th(n, eom, t, y, stepsize)
     if(mod(i, 100).eq.0)   write(*, "(7G14.5)") t*1.e6, y(1:3)*1.e3, y(4:6) !!output unit: microsecond, mm, m/s !test
  enddo
  
end program pt


