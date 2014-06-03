module basic_utils
  implicit none
  integer,parameter::IB=kind(1)
  integer(IB),parameter::dl=kind(1.d0)
  integer(IB),parameter::dlc=kind((1.d0,1.d0))
  integer(IB),parameter::sp=kind(1.)
  Character,Parameter::const_backslash = Char(92)
  character,parameter::const_slash = Char(47)
  Character,Parameter::const_backspace = Char(8)
  Character,Parameter::const_tab = Char(9)
  Character,Parameter::const_newline = Char(10)
  Character,Parameter::const_vertical_tab = Char(11)
  Character,Parameter::const_newpage = Char(12)
  Character,Parameter::const_carriage_return = Char(13)
  real(dl),parameter:: const_pi = 3.14159265358979323846264338327950288d0
  real(dl),parameter:: const_pio2 = const_pi/2._dl
  real(dl),parameter:: const_pio4 = const_pi/4._dl
  real(dl),parameter:: const_ln2 = 0.6931471805599453094172321d0 
  real(dl),parameter:: const_ln10 = 2.302585092994045684017991d0
  real(dl),PARAMETER:: const_LnPi=1.144729885849400174143427351d0
  real(dl),PARAMETER:: const_LogPi=0.4971498726941338543512682883d0
  real(dl),parameter:: const_sqrt2 = 1.4142135623730950488016887d0
  real(dl),PARAMETER:: const_sqrt3 = 1.73205080756887729352744634d0
  real(dl),parameter:: const_sqrt5 = 2.236067977499789696409174d0
  real(dl),parameter:: const_sqrt6 = const_sqrt2*const_sqrt3
  real(dl),parameter:: const_sqrt7 = 2.645751311064590590502d0
  real(dl),parameter:: const_sqrt8 = 2.d0*const_sqrt2
  real(dl),parameter:: const_sqrt10 = const_sqrt2 * const_sqrt5
  real(dl),parameter:: const_sqrt11 = 3.316624790355399849115d0
  real(dl),parameter:: const_sqrt12 = 2.d0*const_sqrt3
  real(dl),parameter:: const_sqrt13 = 3.605551275463989293119d0
  real(dl),parameter:: const_sqrt14 = const_sqrt2*const_sqrt7
  real(dl),parameter:: const_sqrt15 = const_sqrt3*const_sqrt5
  real(dl),parameter:: const_8pi = const_pi*8._dl
  real(dl),parameter:: const_4pi = const_pi*4._dl
  real(dl),parameter:: const_4piby3 = const_4pi/3._dl
  real(dl),parameter:: const_3by4pi = 3._dl/const_4pi
  real(dl),parameter:: const_2pi = const_pi*2._dl
  real(dl),parameter:: const_pi2 = const_pi ** 2
  real(dl),parameter:: const_pi4 = const_pi2 ** 2
  real(dl),parameter:: const_2pi2 = const_pi2 * 2._dl
  real(dl),parameter:: const_8pi3 = const_2pi ** 3
  real(dl),parameter:: const_7pi4by120 = const_pi4 * (7.d0/120.d0)
  real(dl),parameter:: const_third = 1._dl/3._dl
  real(dl),parameter:: const_two_thirds = 2._dl/3._dl
  real(dl),parameter:: const_four_thirds = 4._dl/3._dl
  real(dl),parameter:: const_sqrtpi = 1.7724538509055160272981674833411_dl
  real(dl),parameter:: const_EulerC=0.57721566490153286060651209_dl
  real(dl),parameter:: const_Riemannzeta3 = 1.2020569031595942853997_dl
  real(dl),parameter:: const_Riemannzeta5  = 1.0369277551433699263313_dl
  real(dl),parameter:: const_Riemannzeta7  = 1.0083492773819228268397_dl
  real(dl),parameter:: const_Riemannzeta9 = 1.00200839282608221441785_dl
  real(dl),PARAMETER:: const_fullsky_degrees = 41252.96125   !!4pi/degree^2
  real(dl),parameter:: const_degree = const_pi/180.d0
  real(dl),parameter:: const_arcmin = const_degree/60.d0
  real(dl),parameter:: const_sigmabyfwhm = 1.d0/sqrt(8.d0*const_ln2)
  real(dl),parameter:: const_arcsec = const_arcmin/60.d0
  real(dl),parameter:: const_chbyMpcH0 = 2997.92458d0
  real(dl),parameter :: const_NewtonG = 6.67428e-11_dl
  real(dl),parameter:: const_planck = 6.626068d-34
  real(dl),parameter:: const_hbar = const_planck/const_2pi
  real(dl),parameter:: const_boltzmann = 1.3806504d-23
  real(dl),parameter:: const_K2GHz = const_boltzmann/(const_planck*1.d9) !!
  real(dl),parameter:: const_Mpc_SI = 3.08568025e22_dl
  real(dl),parameter:: const_Msun_SI = 1.98892e30_dl
  real(dl),parameter:: const_SpeedOfLight_SI = 299792458._dl
  real(dl),parameter:: const_c = const_SpeedOfLight_SI
  real(dl),parameter:: const_PlanckMass_SI = 2.17644e-8_dl
  real(dl),parameter:: const_PlanckLength_SI =1.616252e-35_dl
  real(dl),parameter:: const_PlanckTime_SI = 5.39124e-44_dl
  real(dl),parameter:: const_PlanckEnergy_SI = const_PlanckMass_SI * const_SpeedOfLight_SI ** 2
  real(dl),parameter:: const_PlanckTemperature_SI = const_PlanckEnergy_SI / const_boltzmann
  real(dl),parameter:: const_Yr_SI = 3600._dl*24._dl*365.2422_dl !!time
  real(dl),parameter:: const_Gyr_SI = const_Yr_SI * 1.e9_dl !!time
  real(dl),parameter:: const_Lyr_SI = 9.4605284e15_dl !!length
  real(dl),parameter:: const_eV_SI = 1.60217648740e-19_dl
  real(dl),parameter:: const_eVmass_SI = const_eV_SI / const_SpeedOfLight_SI ** 2
  real(dl),parameter:: const_eVlength_SI = const_SpeedOfLight_SI/ (const_eV_SI /const_hbar)
  real(dl),parameter:: const_MeV_SI = const_eV_SI * 1.e6_dl
  real(dl),parameter:: const_MeVmass_SI = const_MeV_SI / const_SpeedOfLight_SI ** 2
  real(dl),parameter:: const_GeV_SI = const_eV_SI * 1.e9_dl
  real(dl),parameter:: const_GeVmass_SI = const_GeV_SI / const_SpeedOfLight_SI ** 2

  real(dl),parameter:: const_fine_structure = 1._dl/137.035_dl
  real(dl),parameter:: const_atomic_mass_unit_SI = 1.66053878283e-27_dl

  real(dl),parameter:: const_electron_mass_SI = const_atomic_mass_unit_SI / 1822.8884845_dl !!0.51099892811 * const_MeVmass_SI

  real(dl),parameter:: const_proton_mass_SI = 1.00727646681290 * const_atomic_mass_unit_SI

  real(dl),parameter:: const_E_hydrogen = 13.605698 * const_eV_SI 

  real(dl),parameter:: const_Hydrogen_mass_SI = const_electron_mass_SI + const_proton_mass_SI - const_E_Hydrogen/const_SpeedOfLight_SI**2

  real(dl),parameter:: const_massratio_He_H =  3.9715_dl

  real(dl),parameter:: const_rhocritbyh2_SI = 3._dl*(1.e5_dl/const_SpeedOfLight_SI)**2*(const_PlanckMass_SI/const_Mpc_SI**2/const_PlanckLength_SI)/const_8pi !!1.878e-26 !!mass/volume, note that this is NOT energy per volume
  real(dl),parameter:: const_hbyH0_SI = const_Mpc_SI / 1.e5_dl !!Cosmic age unit
  real(dl),parameter:: const_hbyH0_Gyr = const_hbyH0_SI / const_Gyr_SI
  real(dl),parameter:: const_chbyH0_Mpc = const_chbyMpcH0  !!alias
  real(dl),parameter:: const_chbyH0_SI = const_chbyH0_Mpc * const_Mpc_SI  !!alias
  real(dl),parameter:: const_H0byh_SI = 1./const_hbyH0_SI
  real(dl),parameter:: const_H0byh_eV = const_hbar/const_eV_SI * const_H0byh_SI

  real(dl),parameter:: const_blackbody_alpha_SI = (const_pi2/30._dl)*(const_boltzmann) / (const_PlanckTemperature_SI * const_PlanckLength_SI)**3  !!blackbody radiation density = alpha * (total g) * (Temperature in Kelvin) ^ 4  [g=1 for each spin degree of boson, g=7/8 for each spin degree of fermion; for photon g_total = 2; for neutrinos g_total = 7/8 * number of species] (result is in SI unit J/m^3)
  real(dl),parameter:: const_Stefan_Boltzmann = 5.670400e-8_dl !!W/m^2/K^4

  real(dl),parameter:: const_sigma_thomson_SI = 6.6524616e-29_dl

  real(dl), parameter :: const_arad_SI = (const_pi**5 * 8./15. ) * const_boltzmann * (const_boltzmann/(const_SpeedOfLight_SI * const_Planck)) **3  
    !7.565914e-16_dl !radiation constant for u=aT^4
  !! = 2.*const_blackboday_alpha_SI (since photon has two spin dof)

  real(dl), parameter :: const_ComptonCT = (8.d0/3.d0) * const_sigma_thomson_SI/(const_electron_mass_SI * const_SpeedOfLight_SI**2) * const_arad_SI  * const_chbyH0_SI  !! ch/H_0 * 8/3 alpha * sigma_T / (m_e c^2)
  real(dl),parameter::const_barssc0 = const_boltzmann / const_hydrogen_mass_SI / const_SpeedOfLight_SI ** 2


  real(dl),parameter:: const_rhocritbyh2_MsunbyMpc3 = 3._dl*(const_PlanckMass_SI/const_Msun_SI)*(const_Mpc_SI/const_PlanckLength_SI) * (1.e5_dl/const_SpeedOfLight_SI)**2 / const_8pi !!2.77467e11
  real(dl),parameter::const_deltac_EdS = (3._dl/5._dl)*(3._dl/2._dl*const_pi)**(2._dl/3._dl) !!spherical collapse critical density contrast

  integer,parameter::factorial01=1
  integer,parameter::factorial02=2
  integer,parameter::factorial03=6
  integer,parameter::factorial04=24
  integer,parameter::factorial05=120
  integer,parameter::factorial06=720
  integer,parameter::factorial07=5040
  integer,parameter::factorial08=40320
  integer,parameter::factorial09=362880
  integer,parameter::factorial10=3628800
  integer,parameter::factorial11=39916800
  integer,parameter::factorial12=479001600

  Integer(IB),parameter::tmp_file_unit = 7
  Integer(IB),parameter::ps_file_unit=8


  Integer(IB),parameter::params_pass_dim_f = 32
  Integer(IB),parameter::params_pass_dim_i = 8
  Integer(IB),parameter::params_pass_dim_l = 8

  Type Params_Pass  !!to avoid using a lot of global variables (which are difficult to parallelize), you may want to write subroutines with a small "parameter set"
     real(dl) fp(params_pass_dim_f)
     integer(IB) ip(params_pass_dim_i)
     logical lp(params_pass_dim_l)
  End type Params_Pass


  interface set_uniform
     module procedure set_uniform_s, set_uniform_d
  end interface set_uniform


  interface findgen
     module procedure findgen_s, findgen_d
  end interface findgen


  Interface Swap
     module procedure swap_float,swap_int, swap_int_array,swap_float_array
  end Interface

  interface is_triangle
     module procedure is_triangle_real, is_triangle_int
  end interface is_triangle

  Interface signln
     module procedure signln_s, signln_v, signln_sp, signln_vp
  end Interface signln

  Interface signexp
     module procedure signexp_s, signexp_v, signexp_sp, signexp_vp
  end Interface signexp


  Interface get_signln
     module procedure get_signln_s, get_signln_v,get_signln_sp, get_signln_vp
  end Interface get_signln

  Interface get_signexp
     module procedure  get_signexp_s, get_signexp_v, get_signexp_sp, get_signexp_vp
  end Interface get_signexp

  logical::global_feedback = .true.
  logical::global_debug_mode = .false.
  integer::global_counter = 0
  logical::global_error = .false.

contains


  subroutine swap_int(x,y)
    integer(IB) x,y,tmp
    tmp=x
    x=y
    y=tmp
  end subroutine swap_int

  subroutine swap_int_array(ix,iy)
    integer(IB),dimension(:),intent(INOUT)::ix
    integer(IB),dimension(:),intent(INOUT)::iy
    integer(IB) tmp(size(ix)) 
    if(size(ix).ne.size(iy)) stop "can not swap arrays with different sizes"
    tmp=ix
    ix=iy
    iy=tmp
  end subroutine swap_int_array

  subroutine swap_float(x,y)
    real(dl) x,y,tmp
    tmp=x
    x=y
    y=tmp
  end subroutine swap_float

  subroutine swap_float_array(xa,xb)
    real(dl),dimension(:),intent(INOUT)::xa,xb
    real(dl) tmp(size(xa))
    if(size(xb).ne.size(xa)) stop "can not swap arrays with different sizes"
    tmp=xa
    xa=xb
    xb=tmp
  end subroutine swap_float_array


  subroutine set_uniform_s(n, x, xstart, xend)
    integer n, i
    real x(n), xstart, xend
    real dx, sumx
    if(n.le.1)then
       if(n.eq.1) x = xstart
       return
    endif
    dx = (xend-xstart)/(n-1._dl)
    x(1)=xstart
    x(n)=xend
    sumx = xstart + xend
    !$omp parallel do
    do i = 2, (n+1)/2
       x(i) = xstart + (i-1) * dx
       x(n+1-i) = sumx - x(i)
    enddo
    !$omp end parallel do
  end subroutine set_uniform_s

  subroutine set_uniform_d(n, x, xstart, xend)
    integer n, i
    real(dl) x(n), xstart, xend
    real(dl) dx, sumx
    if(n.le.1)then
       if(n.eq.1) x = xstart
       return
    endif
    dx = (xend-xstart)/(n-1._dl)
    x(1)=xstart
    x(n)=xend
    sumx = xstart + xend
    !$omp parallel do
    do i = 2, (n+1)/2
       x(i) = xstart + (i-1) * dx
       x(n+1-i) = sumx - x(i)
    enddo
    !$omp end parallel do
  end subroutine set_uniform_d

  subroutine findgen_d(x, xstart,xend)
    real(dl),DIMENSION(:),INTENT(out)::x
    real(dl) xstart,xend
    call set_uniform_d(size(x), x, xstart, xend)
  end subroutine FIndGen_d

subroutine findgen_s(x, xstart,xend)
    real,DIMENSION(:),INTENT(out)::x
    real xstart,xend
    call set_uniform_s(size(x), x, xstart, xend)
  end subroutine FIndGen_s


  subroutine FIndGen_center(x,xstart,xend)
    real(dl),DIMENSION(:),INTENT(out)::x
    real(dl) xstart,xend,dx
    integer(IB) i,n
    n=size(x)
    dx = (xend-xstart)/real(n, dl)
    !$omp parallel do
    do i = 1, n
       x(i) = xstart + (i-0.5_dl) * dx
    enddo
    !$omp end parallel do
  end subroutine FIndGen_center

  subroutine FIndGen_noedge(x,xstart,xend)
    real(dl),DIMENSION(:),INTENT(out)::x
    real(dl) xstart,xend,dx
    integer(IB) i,n
    n=size(x)
    dx = (xend-xstart)/(n+1._dl)
    !$omp parallel do
    do i = 1, n
       x(i) = xstart + i * dx
    enddo
    !$omp end parallel do
  end subroutine FIndGen_noedge

  subroutine FIndGenLn(x,xstart,xend)
    real(dl),DIMENSION(:),INTENT(out)::x
    real(dl) xstart,xend,dx, ln1, ln2
    integer(IB) i,n
    if(xstart.le.0 .or. xend .le.0) call ReturnError("FindGenLn","Cannot do uniform log for negative starting/end", "stop")
    n=size(x)
    ln1 = log(xstart)
    ln2 = log(xend)
    dx = (ln2-ln1)/(n-1._dl)
    x(1)=xstart
    x(n)=xend
    do i = 2, n-1
       x(i) = exp(ln1 + (i-1) * dx)
    enddo
  end subroutine FIndGenLn

  subroutine Debug_Pause(modn)
    integer,save::i = -1
    character(LEN=32) str
    integer,optional::modn
    !$omp critical
    i=i+1
    !$omp end critical
    if(present(modn))then
       if(mod(i,modn).ne.0) return
    endif
    write(*,'(A, I5)') "Program pause point: ", i
100 write(*,'(A)') "Do you want to continue? (Y/N)"
    read(*,*) str
    str=trim(adjustl(str))
    select case(str(1:1))
    case("y", "Y")
       return
    case("n", "N")
       stop "Program stopped as you entered N"
    case default
       goto 100
    end select
  end subroutine Debug_Pause


  Subroutine ReturnError(name, message, action)
    CHARACTER(LEN=*),optional::name
    Character(Len=*),optional::Message
    Character(Len=*),optional::action
    Character::ans
    if(present(name)) WRITE(*,*) "Error/Warning in "//TRIM(name)
    if(present(Message)) write(*,*) "** "//Trim(Message)//" **"
    if(present(ACTION))THEN
       select case(action)
       case("stop","STOP","Stop", "abort", "Abort", "ABORT")
          stop "program terminated"
       case("return","RETURN","Return", "pass", "Pass", "PASS")
          return
       case("wait","Wait","WAIT", "ask", "ASK", "Ask")
100       write(*,*) "Ignore this error and continue? (Y/N)"
          read(*,*) ans
          if(ans == "y" .or. ans == "Y")then
             return
          else
             goto 100
          endif
       end select
    else
       stop
    endif
  end subroutine ReturnError

  subroutine print_bar()
    write(*,*) "***********************************************************"
  end subroutine print_bar

  function systime_sec(reset) result(sec)
    real(dl) sec
    logical,optional::Reset
    real(dl),save::pretime=0._dl 
    integer nowtime, countrate
    call system_clock(nowtime, countrate)
    if(present(Reset))then
       if(Reset)then
          Pretime = NowTime/dble(CountRate)
          sec = 0._dl
          return
       end if
    end if
    sec = NowTime/dble(CountRate)-PreTime
    return 
  end function systime_sec

  subroutine PrtSysTime(Reset)
    logical,optional::Reset
    real(dl),save::pretime=0.d0
    integer(IB) Nowtime,Countrate
    call system_clock(NOWTIME,COUNTRATE)
    If(Present(Reset))then
       Pretime=Nowtime/Real(Countrate)
       If(Reset)Then
          print*,"===== time label reset to zero. ===="
       endif
    else
       Write(*,'(A18,F14.4,A10)') "===== time label: ",Nowtime/Real(Countrate)-pretime," sec ====="
    endif
  end subroutine PrtSysTime

  subroutine init_random()
    integer(IB) nowtime, countrate, s(3), si
    call System_clock(NOWTIME,COUNTRATE)
    s(1) = nowtime
    s(2) = countrate
    s(3) = mod(nowtime, 17)+1
    call random_seed(SIZE = si)
    call random_seed(PUT = s(1:si))
  end subroutine init_random


!!$#ifdef __GFORTRAN__
!!$  function iargc ()
!!$    integer(IB) iargc
!!$    iargc=command_argument_count()
!!$  end function iargc
!!$  
!!$  subroutine getarg(num, res)
!!$    integer(IB), intent(in) :: num
!!$    character(len=*), intent(out) :: res
!!$    integer(IB) l, err
!!$    call get_command_argument(num,res,l,err)
!!$  end subroutine getarg
!!$#endif

  function GetParam(i)
    character(LEN=512) GetParam
    integer(IB), intent(in) :: i
    if (iargc() < i) then
       GetParam = ''
    else
       call getarg(i,GetParam)
    end if
  end function GetParam



  !!%%%%%%%%%%%%%%%%%% GETDIM %%%%%%%%%%%%%%%%%%%%%
  function GETDIM(SUBRTNAME,N1,N2,N3,N4,N5,N6)
    CHARACTER(LEN=*) SUBRTNAME
    integer(IB),INTENT(IN)::N1,N2
    integer(IB),OPTIONAL::N3,N4,N5,N6
    integer(IB) GETDIM
    GETDIM=N1
    if(N1.NE.N2)THEN
       write(*,"(A,2I8)") "Dimension Error:",N1,N2
       CALL ReturnError(subRtName)
       return
    ELSE
       if(present(N3))THEN
          if(N1.NE.N3)THEN
             write(*,"(A,3I8)") "Dimension Error:",N1,N2, N3
             CALL ReturnError(SUBRTNAME)
             return
          Endif
          if(present(N4))THEN
             if(N1.NE.N4)then
                write(*,"(A,4I8)") "Dimension Error:",N1,N2, N3, N4
                call ReturnError(SubRtName)
                return
             endif
             if(Present(N5))then
                if(N5.ne.N1)then
                   write(*,"(A,5I8)") "Dimension Error:",N1,N2, N3, N4, N5
                   call ReturnError(SubRtName)
                   return
                endif
                if(present(N6))then
                   if(N1.ne.N6)then
                      write(*,"(A,6I8)") "Dimension Error:",N1,N2, N3, N4, N5, N6
                      call ReturnError(SubRtName)
                      return
                   endif
                endif
             endif
          endif
       endif
    endif
  end function GETDIM

  function imaxloc(arr)
    real(dl), DIMENSION(:), INTENT(IN) :: arr
    integer(IB) imaxloc
    integer(IB),DIMENSION(1)::imax
    imax=maxloc(arr(:))
    imaxloc=imax(1)
  end function imaxloc

  function iminloc(arr)
    real(dl), DIMENSION(:), INTENT(IN) :: arr
    integer(IB) iminloc
    integer(IB),DIMENSION(1)::imin
    imin=minloc(arr(:))
    iminloc=imin(1)
  end function iminloc


  function leftloc(arr, x)
    integer leftloc
    real(dl) x
    real(dl),dimension(:),intent(IN)::arr
    integer n, rightloc, midloc
    n =size(arr)
    if((arr(1).lt. x .and. arr(n) .lt. x .and. arr(n).ge.arr(1)).or.(arr(1).gt.x .and. arr(n).gt.x .and. arr(1).ge.arr(n)))then
       leftloc = n
       return
    elseif((arr(1).lt. x .and. arr(n) .lt. x .and. arr(n).le.arr(1)).or.(arr(1).gt.x .and. arr(n).gt.x .and. arr(1).le.arr(n)))then
       leftloc = 0
       return
    endif
    leftloc = 1
    rightloc = n
    if(arr(1) .lt. arr(n))then
       do while(rightloc - leftloc .gt. 1)
          midloc = (leftloc+rightloc)/2
          if(arr(midloc).lt.x)then
             leftloc = midloc
          else
             rightloc =midloc
          endif
       enddo
    else
       do while(rightloc - leftloc .gt. 1)
          midloc = (leftloc+rightloc)/2
          if(arr(midloc).gt.x)then
             leftloc = midloc
          else
             rightloc =midloc
          endif
       enddo
    endif
  end function leftloc

  function nearest_index(arr, x) result(ind)
    real(dl) x
    real(dl),dimension(:),intent(in)::arr
    integer ind
    ind = leftloc(arr, x)
    if(abs(arr(ind) - x) .gt. abs(arr(ind+1)-x))then
       ind = ind+1
    endif
  end function nearest_index

  Function OuterProd(A,B)
    real(dl), DIMENSION(:), INTENT(IN) :: a,b
    real(dl), DIMENSION(size(a),size(b)) :: outerprod
    Outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  end function OuterProd

  function lnabs(x)
    real(dl) x, lnabs
    real(dl),parameter::log_zero = -1.e30_dl
    if(x.eq.0.d0)then
       lnabs = log_zero
    else
       lnabs = dlog(dabs(x))
    endif
  end function lnabs

  subroutine multiply_exp(A, x) !!replace A with A*exp(x), A can be negative
    real(dl) A, x
    A = sign(dexp(x+lnabs(A)),A)
  end subroutine multiply_exp

  function is_triangle_real(k1, k2, k3)
    real(dl) k1, k2, k3
    logical is_triangle_real
    is_triangle_real = (k1+k2 .gt. k3 .and. k2+k3 .gt. k1 .and. k1+k3 .gt. k2)
  end function is_triangle_real

  function is_triangle_int(l1, l2, l3)
    integer l1, l2, l3
    logical is_triangle_int
    is_triangle_int = (l1+l2 .ge. l3 .and. l2+l3.ge.l1 .and. l1+l3.ge.l2)
  end function is_triangle_int

  function minus1power(l)
    integer l
    real(dl) minus1power
    if(mod(l,2).eq.0)then
       minus1power = 1.d0
    else
       minus1power = -1.d0
    endif
  end function minus1power
  
  subroutine feedback_print(str,i)
    character(LEN=*),optional:: str
    integer,optional::i
    character(LEN=32) istr
    if(global_feedback)then
       if(present(str))then
          if(present(i))then
             write(istr,*) i
             write(*,*) str//trim(istr)
          else
             write(*,*) str
          endif
       else
          if(present(i))then
             write(istr,*) i
             write(*,*) trim(istr)
          endif
       endif
    endif
  end subroutine feedback_print

  function signln_s(x) result(s)
    real(dl) x, s
    s = sign(dlog(1.d0+dabs(x)), x)
  end function signln_s

  function signln_sp(x) result(s)
    real(sp) x, s
    s = sign(log(1.+abs(x)), x)
  end function signln_sp

  function signln_v(x) result(s)
    real(dl) x(:)
    real(dl) s(size(x))
    s = sign(dlog(1.d0+dabs(x)), x)
  end function signln_v

  function signln_vp(x) result(s)
    real(sp) x(:)
    real(sp) s(size(x))
    s = sign(log(1.+abs(x)), x)
  end function signln_vp


  function signexp_s(x) result(s)
    real(dl) x, s
    s = sign(dexp(dabs(x))-1.d0, x)
  end function signexp_s

  function signexp_v(x) result(s)
    real(dl) x(:)
    real(dl) s(size(x))
    s = sign(dexp(dabs(x))-1.d0, x)
  end function signexp_v

  function signexp_sp(x) result(s)
    real(sp) x, s
    s = sign(exp(abs(x))-1., x)
  end function signexp_sp

  function signexp_vp(x) result(s)
    real(sp) x(:)
    real(sp) s(size(x))
    s = sign(exp(abs(x))-1., x)
  end function signexp_vp

  subroutine  get_signln_s(x)
    real(dl),intent(INOUT):: x
    x = sign(dlog(1.d0+dabs(x)), x)
  end subroutine get_signln_s

  subroutine  get_signln_v(x)
    real(dl),dimension(:),intent(INOUT)::x
    x = sign(dlog(1.d0+dabs(x)), x)
  end subroutine get_signln_v

  subroutine  get_signln_sp(x)
    real(sp),intent(INOUT):: x
    x = sign(log(1.+abs(x)), x)
  end subroutine get_signln_sp

  subroutine  get_signln_vp(x)
    real(sp),dimension(:),intent(INOUT)::x
    x = sign(log(1.+abs(x)), x)
  end subroutine get_signln_vp

  subroutine get_signexp_s(x)
    real(dl) x
    x = sign(dexp(dabs(x))-1.d0, x)
  end subroutine get_signexp_s

  subroutine get_signexp_v(x)
    real(dl),dimension(:),intent(INOUT)::x
    x = sign(dexp(dabs(x))-1.d0, x)
  end subroutine get_signexp_v

  subroutine get_signexp_sp(x)
    real(sp) x
    x = sign(exp(abs(x))-1., x)
  end subroutine get_signexp_sp

  subroutine get_signexp_vp(x)
    real(sp),dimension(:),intent(INOUT)::x
    x = sign(exp(abs(x))-1., x)
  end subroutine get_signexp_vp

  function query_yn(question) result(ans)
    logical ans
    character(LEN=*) question
    character(LEN=1) input
    write(*,"(A)") trim(question)//" (enter Y or N)"
    read(*,*) input
    select case(input)
    case("Y","y")
       ans = .true.
    case default
       ans = .false.
    end select
  end function query_yn
  
  function true1_false0(istrue)
    logical istrue
    real(dl) true1_false0
    if(istrue)then
       true1_false0 = 1.d0
    else
       true1_false0 = 0.d0
    end if
  end function true1_false0

  subroutine get_polyvalue(p, n, x, y)
    integer n
    real(dl) p(n)
    real(dl) x, y
    integer j
    y = p(n)
    do j= n-1,1,-1
       y= y*x + p(j)
    enddo
  end subroutine get_polyvalue

  function simple_polyvalue(p, n, x) result(y)
    integer n
    real(dl) p(n)
    real(dl) x, y
    integer j
    y = p(n)
    do j= n-1,1,-1
       y= y*x + p(j)
    enddo    
  end function simple_polyvalue


  subroutine vector_cross_product(v1, v2, v3)
    real(dl),dimension(3)::v1,v2,v3
    v3 = (/ v1(2)*v2(3) - v1(3)*v2(2),  v1(3)*v2(1)-v1(1)*v2(3), v1(1)*v2(2)-v1(2)*v2(1) /)
  end subroutine vector_cross_product

end module basic_utils
