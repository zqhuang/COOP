
  function coop_cosmology_firstorder_late_damp_factor(this, k, tau) result(Gam)
    class(coop_cosmology_firstorder)::this
    COOP_REAL::k, Gam, tau
    if(coop_firstorder_optimize)then
       Gam = (1.d0 - tanh((k*tau-max(1000.d0, k*this%tau_late))/50.d0))/2.d0
    else
       Gam = 1.d0
    endif
  end function coop_cosmology_firstorder_late_damp_factor

  subroutine coop_cosmology_firstorder_init_source(this, m)
    class(coop_cosmology_firstorder)::this
    COOP_INT :: m
    this%source(m)%distlss = this%distlss
    this%source(m)%m = m
    this%source(m)%nsrc = coop_num_sources(m)
    this%source(m)%nsaux = coop_num_saux(m)
    call this%set_source_tau(this%source(m), coop_source_tau_step_factor(m))
    call this%set_source_k(this%source(m), coop_source_k_n(m), coop_source_k_weight(m))

    if(allocated(this%source(m)%s))then
       if(size(this%source(m)%s, 1) .ne. this%source(m)%ntau .or. size(this%source(m)%s, 2) .ne. this%source(m)%nk)then
          deallocate(this%source(m)%s, this%source(m)%s2, this%source(m)%saux)
          allocate(this%source(m)%s(this%source(m)%nsrc, this%source(m)%nk , this%source(m)%ntau), this%source(m)%s2(this%source(m)%nsrc, this%source(m)%nk, this%source(m)%ntau), this%source(m)%saux(this%source(m)%nsaux, this%source(m)%nk, this%source(m)%ntau) )
       endif
    else
       allocate(this%source(m)%s(this%source(m)%nsrc, this%source(m)%nk , this%source(m)%ntau), this%source(m)%s2(this%source(m)%nsrc, this%source(m)%nk, this%source(m)%ntau), this%source(m)%saux(this%source(m)%nsaux, this%source(m)%nk, this%source(m)%ntau) )
    endif
  end subroutine coop_cosmology_firstorder_init_source

  function coop_cosmology_firstorder_Clzetazeta_at_R(this, l, r) result(Cl)
    class(coop_cosmology_firstorder)::this
    !!this computeys 4 \pi \int_0^\infty |j_l(kr)|^2 (k^3P(k)/(2\pi^2)) d\ln k
    COOP_REAL::Cl, r
    COOP_INT l
    Cl = (coop_pi**2/2.d0)* this%psofk(1.d0/r) * 2.d0 ** ( this%ns - 1.d0) * exp(log_gamma(3.d0 - this%ns )+log_gamma(l + (this%ns - 1.d0)/2.d0)-log_gamma(l+(5.d0-this%ns)/2.d0)-2.d0*log_gamma((4.d0-this%ns)/2.d0))
  end function Coop_cosmology_firstorder_Clzetazeta_at_R


  function coop_cosmology_firstorder_Clzetazeta(this, l, r1, r2) result(Cl)
    !!this computes 4 \pi \int_0^\infty j_l(kr_1) j_l(kr_2) (k^3P(k)/(2\pi^2)) d\ln k
    class(coop_cosmology_firstorder)::this
    COOP_REAL Cl, r1
    COOP_REAL, optional::r2
    COOP_INT l
    if(.not. present(r2))then
       Cl = this%Clzetazeta_at_R(l, r1)
    else
       if(r1.eq.r2)then
          Cl = this%Clzetazeta_at_R(l, r1)
       else
          cl = this%psofk(2.d0/(r1+r2))*coop_4pi*coop_sphericalBesselCross(l,l,r1,r2,this%ns)
       endif
    endif
  end function Coop_cosmology_firstorder_Clzetazeta

  subroutine coop_cosmology_firstorder_source_intbypart(source, ik)
    COOP_REAL,parameter::k_trunc = 2.d0  !!for tiny k the integrated-by-part term is negligible.
    class(coop_cosmology_firstorder_source)::source
    COOP_REAL s2(source%ntau), sum1
    COOP_INT ik, i
    if(source%k(ik).lt. k_trunc/5.d0 )return
    source%saux(1, ik, 1) = 0.d0
    source%saux(1, ik, source%ntau) = 0.d0
    s2(1) = ((source%saux(1, ik, 1) -  source%saux(1, ik, 2))/(source%dtau(1)+source%dtau(2)) +  (source%saux(1, ik, 1) )/(source%dtau(1)*2.d0))*(2.d0/source%dtau(1))
    s2(source%ntau) = (source%saux(1, ik, source%ntau)/(source%dtau(source%ntau)*2.d0) +  (source%saux(1, ik, source%ntau) -  source%saux(1, ik, source%ntau-1))/(source%dtau(source%ntau)+source%dtau(source%ntau-1)))*(2.d0/source%dtau(source%ntau))
    do i = 2, source%ntau-1
       s2(i) = ((source%saux(1, ik, i) -  source%saux(1, ik, i+1))/(source%dtau(i)+source%dtau(i+1)) +  (source%saux(1, ik, i) -  source%saux(1, ik, i-1))/(source%dtau(i)+source%dtau(i-1)))*(2.d0/source%dtau(i))
    enddo
    source%saux(1, ik, :) = s2*(-(tanh(source%k(ik)/k_trunc))**8)/source%k(ik)**2
    select case(source%m)
    case(0)
       source%s(1, ik, :) = source%s(1, ik, :) +  source%saux(1, ik, :)
    case(1)
       call coop_tbw("vector intbypart")
    case(2)
       source%s(2, ik, :) = source%s(2, ik, :) +  source%saux(1, ik, :)
    case default
       call coop_return_error("source_intbypart", "unknown m", "stop")
    end select
  end subroutine coop_cosmology_firstorder_source_intbypart



!!============================================================================

  subroutine coop_cosmology_firstorder_set_klms(this)
    class(coop_cosmology_firstorder)::this
    COOP_INT k, l, m, s
    if(this%klms_done) return
    do s = 0, coop_pert_default_smax
       do m = 0,coop_pert_default_mmax
          do l = 0, coop_pert_default_lmax
             if(m.ge.l .or. s.ge.l)then
                this%klms(l,m,s) = 0.d0
             else
                this%klms(l, m, s) = sqrt(dble(l**2-m**2)*dble(l**2-s**2))/l
             endif
             this%klms_by_2lp1(l, m, s) = this%klms(l, m, s)/(2*l+1.d0)
             this%klms_by_2lm1(l, m, s) = this%klms(l, m, s)/(2*l-1.d0)
          enddo
       enddo
    enddo
    do l=1, coop_pert_default_lmax
       this%fourbyllp1(l) = 4.d0/l/(l+1.d0)
    enddo
    this%klms_done = .true.
  end subroutine coop_cosmology_firstorder_set_klms


  function coop_cosmology_firstorder_cs2bofa(this, a) result(cs2bofa)
    class(coop_cosmology_firstorder)::this
    COOP_REAL a, cs2bofa
    cs2bofa = this%species(this%index_baryon)%fcs2%eval(a)
  end function coop_cosmology_firstorder_cs2bofa

  function coop_cosmology_firstorder_Tbofa(this, a) result(Tbofa)
    class(coop_cosmology_firstorder)::this
    COOP_REAL a, Tbofa
    Tbofa = this%Tb%eval(a)
  end function coop_cosmology_firstorder_Tbofa


  subroutine coop_cosmology_firstorder_source_free(this)
    class(coop_cosmology_firstorder_source)::this
    if(allocated(this%k))deallocate(this%k, this%dk, this%index_tc_off, this%kop)
    this%nk = 0
    if(allocated(this%tau))deallocate(this%tau, this%chi, this%dtau, this%a, this%tauc, this%lna, this%omega_rad, this%vis, this%omega_de)
    this%ntau = 0
    if(allocated(this%s))deallocate(this%s, this%s2, this%saux)
    if(allocated(this%k_dense))deallocate(this%k_dense, this%ws_dense, this%wt_dense, this%ps_dense)
    COOP_DEALLOC(this%Cls)
    COOP_DEALLOC(this%Cls_lensed)
    call this%trans%free()
  end subroutine coop_cosmology_firstorder_source_free

  subroutine coop_cosmology_firstorder_free(this)
    class(coop_cosmology_firstorder)::this
    COOP_INT m, i
    this%index_baryon = 0
    this%index_cdm = 0
    this%index_de = 0
    this%index_wdm = 0
    this%index_radiation = 0
    this%index_nu = 0
    this%index_massivenu = 0
    do m=0,2
       call this%source(m)%free()
    enddo
    call this%Ps%free()
    call this%Pt%free()
    call this%pssmooth%free()
    call this%ptsmooth%free()
    call this%Xe%free()
    call this%eKappa%free()
    call this%vis%free()
    call this%Tb%free()
    call this%fdis%free
    call this%ftime%free
    call this%faoftau%free
    do i= 1, this%num_species
       call this%species(i)%free
    enddo
    this%num_species = 0
    this%omega_k_value = 1.
    this%need_setup_background = .true.
#if DO_EFT_DE
       call this%f_Mpsq%free()
       call this%f_alpha_M%free()       
       call this%f_alpha_H%free()
       call this%f_alpha_T%free()
       call this%f_alpha_B%free()       
       call this%f_alpha_K%free()
#endif 
  end subroutine coop_cosmology_firstorder_free


  subroutine coop_cosmology_firstorder_set_power(this, power, args)
!!the format of primordial power
!!subroutine power(k/k_pivot, ps, pt, cosmology, args)
!!input kMpc and args, output ps and pt    
    integer,parameter::n = 100000
    class(coop_cosmology_firstorder)::this
    external power
    interface
       subroutine power(krat, Pscalar, Ptensor, arguments)
	 use coop_wrapper_background
	 use coop_pertobj_mod
         use coop_cl_indices_mod
         use coop_lensing_mod
         COOP_REAL::krat, Pscalar, Ptensor
         type(coop_arguments)::arguments
       end subroutine power
    end interface
    type(coop_arguments) args
    COOP_REAL k(n), ps(n), pt(n), dlnk
    COOP_INT i, nsteps
    this%k_pivot = this%kMpc_pivot/this%H0Mpc()
    call coop_set_uniform(n, k, coop_power_kmin, coop_power_kmax, logscale = .true.)
    dlnk = log(k(2)/k(1))
    !$omp parallel do
    do i=1, n
       call power(k(i)/this%k_pivot, ps(i), pt(i), args)
    end do
    !$omp end parallel do
    if(any(ps.lt. 0.d0)) stop "Primordial scalar power spectrum cannot be <= 0"
    if(any(pt .lt. 0.d0)) stop "Primordial tensor power spectrum cannot be < 0"
    this%has_tensor = any(pt .gt. 1.d-15)
    call this%ps%init(n = n, xmin = k(1), xmax = k(n), f = ps, xlog = .true., ylog = all(ps .gt. 0.d0), fleft = ps(1), fright = ps(n), check_boundary = .false., method = COOP_INTERPOLATE_LINEAR, name = "ScalarPower")
    if(this%has_tensor)then
       call this%pt%init(n = n, xmin = k(1), xmax = k(n), f = pt, xlog = .true., ylog = all(pt .gt. 0.d0), fleft = pt(1), fright = pt(n), check_boundary = .false., method = COOP_INTERPOLATE_LINEAR, name = "TensorPower")
    else
       call this%pt%init_polynomial((/ 0.d0 /), name = "TensorPower")
    endif
    this%As = this%psofk(this%k_pivot)
    if(this%ps%ylog)then
       this%ns = this%ps%derivative_bare(log(this%k_pivot)) + 1.d0
       this%nrun = this%ps%derivative2_bare(log(this%k_pivot))
    else
       this%ns = this%ps%derivative_bare(log(this%k_pivot))/this%As + 1.d0
       this%nrun = this%ps%derivative2_bare(log(this%k_pivot))/this%As - (this%ns-1.d0)**2
    endif

    if(this%has_tensor)then
       this%r =  this%ptofk(this%k_pivot)/this%As
       if(this%r .eq. 0.d0)then
          this%nt = 0.d0
       else
          this%nt = this%pt%derivative(this%k_pivot)*this%k_pivot/(this%r*this%As)
       endif
    else
       this%r = 0.d0
       this%nt = 0.d0
    endif
    !!compute the smoothed power spectrum (for stability of CMB power integrator )
    call coop_set_uniform(n, k, coop_power_kmin, coop_power_kmax)
    !$omp parallel do
    do i=1, n
       call power(k(i)/this%k_pivot, ps(i), pt(i), args)
    end do
    !$omp end parallel do
    nsteps = nint(1.d0/this%distlss/(k(2)-k(1))) !!smooth out any features below delta_l = 1 resolution
    if(nsteps .gt. 0)then
       call coop_smooth_data(n, ps, nsteps, logscale = all(ps.gt.0.d0) )
       if(this%has_tensor)then
          call coop_smooth_data(n, pt, nsteps, logscale = all(pt.gt.0.d0))
       endif
    endif
    call this%pssmooth%init(n = n, xmin = k(1), xmax = k(n), f = ps, ylog = all(ps.gt.0.d0), fleft = ps(1), fright = ps(n), check_boundary = .false., method = COOP_INTERPOLATE_LINEAR, name = "SmoothedScalarPower")
    if(this%has_tensor)then
       call this%ptsmooth%init(n = n, xmin = k(1), xmax = k(n), f = pt, ylog=all(pt.gt.0.d0), fleft = pt(1), fright = pt(n), check_boundary = .false., method = COOP_INTERPOLATE_LINEAR, name = "SmoothedTensorPower")
    else
       call this%ptsmooth%init_polynomial((/ 0.d0 /), name = "TensorPower")
    endif

  end subroutine coop_cosmology_firstorder_set_power

  subroutine coop_cosmology_firstorder_set_xe(this)
    class(coop_cosmology_firstorder)::this
    integer i, m
    integer, parameter::n = 8192
    COOP_REAL ekappa(n), a(n), dkappadinva(n), dkappadtau(n), vis(n)
    call this%setup_background()
    call this%set_klms()
    this%Omega_b = O0_BARYON(this)%Omega    
    this%ombm2h2 = this%Omega_b * this%h()**2 * this%Mpsq0
    this%Omega_c = O0_CDM(this)%Omega
    this%omcm2h2 = this%Omega_c*O0_CDM(this)%rhoa3_ratio(1.d-10)*this%h()**2* this%Mpsq0
    this%Omega_g = O0_RADIATION(this)%Omega
    this%Rbya = 0.75d0*this%omega_b/this%omega_g
    if(this%index_massiveNu .ne. 0)then
       call coop_fermion_get_lnam(log(O0_MASSIVENU(this)%Omega/O0_MASSIVENU(this)%Omega_massless), this%mnu_by_Tnu)
       this%mnu_by_Tnu = exp(this%mnu_by_Tnu)
       this%Omega_nu = O0_NU(this)%Omega + O0_MASSIVENU(this)%Omega_massless
       this%Omega_massivenu = O0_MASSIVENU(this)%Omega - O0_MASSIVENU(this)%Omega_massless
    else
       this%Omega_nu = O0_NU(this)%Omega 
       this%mnu_by_Tnu = 0.d0
       this%Omega_massivenu = 0.d0
    endif

    this%Omega_m = this%Omega_b + this%Omega_c + this%Omega_massivenu !!
    this%Omega_r =  this%Omega_g + this%Omega_nu 
    this%bbks_keq =  (0.073*coop_SI_c/1.d5) * this%h_value * this%omega_m * exp(- this%omega_b - sqrt(2.*this%h_value)*this%omega_b/this%omega_m)

    this%tau_eq = this%conformal_time(this%a_eq)
    this%tau_late = this%conformal_time(this%a_eq*5.d0)
    
    this%dkappadtau_coef = this%species(this%index_baryon)%Omega * this%h() * coop_SI_sigma_thomson * (coop_SI_rhocritbyh2/coop_SI_c**2) * coop_SI_hbyH0 * coop_SI_c/ coop_SI_m_H * (1.d0 - this%YHe())*this%Mpsq0
    if(this%do_reionization)then
       this%reionFrac = 1.d0 + this%YHe()/(coop_m_He_by_m_H * (1.d0-  this%YHe()))
    else
       this%reionFrac = 0.d0
    endif

    call this%set_zre_from_optre()
    select type(this)
    type is(coop_cosmology_firstorder)
       call coop_recfast_get_xe(this, this%xe, this%Tb, this%reionFrac, this%zre, this%deltaz)
    class default
       stop "For compatibility with lower versions of gfortran set_xe only works with type coop_cosmology_firstorder"
    end select
    call coop_set_uniform(n, a, log(coop_visibility_amin), log(coop_scale_factor_today))

    a = exp(a)
    !$omp parallel do
    do i=1, n
       dkappadtau(i) = this%dkappadtau(a(i))
       dkappadinva(i)= - dkappadtau(i)/this%Hratio(a(i))
    enddo
    !$omp end parallel do
    ekappa(n) = 0.d0
    do i=n-1,1, -1
       ekappa(i) = ekappa(i+1) + (dkappadinva(i+1) + dkappadinva(i))*(0.5d0/a(i) - 0.5d0/a(i+1))
       if(ekappa(i) .lt. -200.d0)then
          ekappa(1:i) = -200.d0
          exit
       endif
    enddo
    ekappa = exp(ekappa)
    vis = dkappadtau*ekappa
    call this%ekappa%init(n = n, xmin = a(1), xmax = a(n), f = ekappa, xlog = .true., ylog = .true., fleft = ekappa(1), fright = ekappa(n), check_boundary = .false.)
    call this%vis%init(n = n, xmin = a(1), xmax = a(n), f = vis, xlog = .true., ylog = .true., fleft = vis(1), fright = vis(n), method = COOP_INTERPOLATE_QUADRATIC, check_boundary = .false.)

    this%arecomb =this%vis%maxloc()
    this%zrecomb = 1.d0/this%arecomb - 1.d0
    this%taurecomb = this%tauofa(this%arecomb)
    this%maxvis = this%vis%eval(1.d0/(1.d0+this%zrecomb))
    this%zrecomb_start = this%zrecomb+30.d0
    this%arecomb_start = 1.d0/(1.d0+this%zrecomb_start)
    call this%get_z_drag(this%z_drag)
    call this%get_z_star(this%z_star)    
    this%r_drag = coop_integrate(dsoundda, 1.d-9, 1.d0/(1.d0+this%z_drag))
    this%r_star = coop_integrate(dsoundda, 1.d-9, 1.d0/(1.d0+this%z_star))    
    
    do while(this%vis%eval(1.d0/(this%zrecomb_start+1.d0))/this%maxvis .gt. 1.d-4 .and. this%zrecomb_start .lt. 2.d4)
       this%zrecomb_start = this%zrecomb_start + 30.d0
    enddo

    this%zrecomb_end = this%zrecomb - 10.d0
    do while(this%vis%eval(1.d0/(this%zrecomb_end+1.d0))/this%maxvis .gt. 1.d-3 .and. this%zrecomb_end .gt. 100.d0)
       this%zrecomb_end = this%zrecomb_end  - 10.d0      
    enddo

    this%distlss = this%comoving_angular_diameter_distance(1.d0/(1.d0+this%zrecomb))
    this%tau0 = this%conformal_time(coop_scale_factor_today)


  contains

    function dsoundda(a) result(dsda)
      COOP_REAL a , dsda
      dsda = this%dsoundda(a)
    end function dsoundda
    
  end subroutine coop_cosmology_firstorder_set_xe

  subroutine coop_cosmology_firstorder_set_user_defined_power(this, upar)
    class(coop_cosmology_firstorder)::this
    COOP_REAL upar(coop_num_user_defined_params)
    type(coop_arguments):: args
    call args%init( r = upar )
    select case(this%pp_genre)
    case(100)
       call this%set_power(coop_user_defined_primordial_power_100, args)
    case(101)
       call this%set_power(coop_user_defined_primordial_power_101, args)
    case(102)
       call this%set_power(coop_user_defined_primordial_power_102, args)
    case(103)
       call this%set_power(coop_user_defined_primordial_power_103, args)
    case(104)
       call this%set_power(coop_user_defined_primordial_power_104, args)
    case(105)
       call this%set_power(coop_user_defined_primordial_power_105, args)
    case(106)
       call this%set_power(coop_user_defined_primordial_power_106, args)
    case(107)
       call this%set_power(coop_user_defined_primordial_power_107, args)
    case(108)
       call this%set_power(coop_user_defined_primordial_power_108, args)
    case(109)
       call this%set_power(coop_user_defined_primordial_power_109, args)
    case(110)
       call this%set_power(coop_user_defined_primordial_power_110, args)
    case default
       write(*,*) "pp_genre = "//COOP_STR_OF(this%pp_genre)
       Write(*,*) "Unknown pp genre. For standard power spectrum use pp_genre = 0; for user defined power spectrum use pp_genre = 101 - 110. Modify include/user_defined_primordial_power.h to define your model."
       stop
    end select
    call args%free()

  end subroutine coop_cosmology_firstorder_set_user_defined_power

  subroutine coop_cosmology_firstorder_set_standard_power(this, As, ns, nrun, r, nt, inflation_consistency)
    class(coop_cosmology_firstorder)::this
    COOP_REAL kpivot, As, ns, nrun, r, nt
    type(coop_arguments):: args
    logical,optional::inflation_consistency
    if(present(inflation_consistency))then
       if(inflation_consistency)then
          call args%init( r = (/ As, ns, nrun, r, -r/8.d0 /) )
       else
          call args%init( r = (/ As, ns, nrun, r, nt /) )
       endif
    else
       call args%init( r = (/ As, ns, nrun, r, nt /) )
    endif
    call this%set_power(coop_cosmology_firstorder_standard_power, args)
    call args%free()
  end subroutine coop_cosmology_firstorder_set_standard_power

  subroutine coop_cosmology_firstorder_standard_power(kbykpiv, ps, pt, args)
    COOP_REAL kbykpiv, ps, pt, lnkbykpiv
    type(coop_arguments) args
    lnkbykpiv = log(kbykpiv)
    ps = args%r(1) * exp((args%r(2)-1.d0 + args%r(3)/2.d0*lnkbykpiv)*lnkbykpiv)
    pt = args%r(4)*args%r(1)*exp(lnkbykpiv*(args%r(5)))
  end subroutine coop_cosmology_firstorder_standard_power

  subroutine coop_cosmology_firstorder_source_k2kop(source, k, kop)
    class(coop_cosmology_firstorder_source)::source
    COOP_REAL::k, kop
    kop = k**coop_source_k_index*source%kweight + log(k)
  end subroutine coop_cosmology_firstorder_source_k2kop

  subroutine coop_cosmology_firstorder_source_kop2k(source, kop, k, dkop, dk)
    class(coop_cosmology_firstorder_source)::source
    COOP_REAL kop, k, kalpha
    COOP_REAL, optional::dkop, dk
    call coop_source_kop2k_noindex(kop*coop_source_k_index, source%kweight*coop_source_k_index, kalpha)
    k = kalpha**(1.d0/coop_source_k_index)
    if(present(dkop) .and. present(dk))then
       dk = dkop * k / (1.d0 + source%kweight * coop_source_k_index * kalpha )
    endif
  end subroutine coop_cosmology_firstorder_source_kop2k

  subroutine coop_source_kop2k_noindex(kop, weight, k)
    !!solve the equation k * weight + log(k) = kop
    COOP_REAL kop, weight, k, kmin, kmax, kmid
    if(kop .gt. weight)then
       kmax = kop/weight
       kmin = (kop+1.d0)/(weight+1.d0)
    else
       kmax = exp(kop)
       kmin = max(kop/weight, 1.d-3)
    endif
    do while((kmax-kmin)/kmin .gt. 1.d-8)
       kmid = (kmax+kmin)/2.d0
       if(kmid*weight + log(kmid) .gt. kop)then
          kmax = kmid
       else
          kmin = kmid
       endif
    enddo
    k = (kmax+kmin)/2.d0
  end subroutine coop_source_kop2k_noindex
  
  function coop_cosmology_firstorder_psofk(this, k) result(ps)
    COOP_REAL k, ps
    class(coop_cosmology_firstorder)::this
    ps = this%ps%eval(k)
  end function coop_cosmology_firstorder_psofk


  function coop_cosmology_firstorder_ptofk(this, k) result(pt)
    COOP_REAL k, pt
    class(coop_cosmology_firstorder)::this
    pt = this%pt%eval(k)
  end function coop_cosmology_firstorder_ptofk


!!slow way to calculate theta. You do not need to initialize background before calling this function.
  function coop_cosmology_firstorder_cosmomc_theta(this) result(theta)
    class(coop_cosmology_firstorder)::this
    COOP_REAL zstar, theta, astar, rs, da, omch2_eff
#if DO_COUPLED_DE
    omch2_eff = O0_CDM(this)%rhoa3_ratio(astar)*this%omcm2h2
    zstar = coop_zrecomb_fitting(this%ombm2h2, omch2_eff)
#else
    zstar = coop_zrecomb_fitting(this%ombm2h2, this%omcm2h2)
#endif    
    astar = 1.d0/(1.d0+zstar)    
    rs = coop_integrate(dsoundda, 1.d-8, astar, 1.d-7)  !!to be consistent with CosmoMC
    da = coop_r_of_chi(coop_integrate(dtauda, astar, 1.d0, 1.d-7), this%Omega_k())
    theta  = rs/da
  contains

    function dsoundda(a) result(dsda)  !!use approximations, to be consistent with CosmoMC
      COOP_REAL a, dsda, R
      R = this%ombm2h2 * 3.d4*a
      dsda = dtauda(a)/sqrt(3.d0*(1.d0+R))
    end function dsoundda

    function dtauda(a) 
      COOP_REAL a, dtauda
      dtauda = 1.d0/this%Hasq(a)
    end function dtauda
  end function coop_cosmology_firstorder_cosmomc_theta
	    

  subroutine coop_cosmology_firstorder_set_source_tau(this, source, step_factor)
    class(coop_cosmology_firstorder_source)::source
    class(coop_cosmology_firstorder)::this
    COOP_REAL step_factor, viscut
    !!basic stepsize
    COOP_REAL,parameter::step_ini = 7.d-4
    COOP_INT,parameter::nmax = 16384
    COOP_REAL, parameter::a_factor = 1.05d0
    COOP_REAL step
    COOP_REAL a(nmax), tau(nmax), aend, tautmp
    COOP_INT i, n, nw, iw, j
    logical::check_input, do_inds
    if(step_factor .lt. 0.01d0) stop "Error in set_source_tau: step_factor must > 0.01"
    n = 1
    aend = 1.d0/(1.d0+this%zrecomb_start)
    a(n) = min(aend/2.d0, max(aend/(3.d0*exp((coop_feedback_level-1)*2.d0)), coop_min_scale_factor))
    tau(n) = this%tauofa(a(n))
    
    step = step_ini

    do while(a(n).lt. aend)
       n =n+1
       if(n.gt.nmax) call coop_return_error("set_source_tau", "nmax<n", "stop")
       a(n) = a(n-1) * a_factor
       tau(n) = tau(n-1) + step*step_factor
       tautmp = this%tauofa(a(n))
       if(tautmp .lt. tau(n))then
          tau(n) = tautmp
       else
          a(n) = this%aoftau(tau(n))
       endif
    enddo

    call expand_array(amax = 1.d0/(1.d0 + max(this%zrecomb+400.d0, 1.d0)), stepfac = 0.98d0, stepmin = step_ini*0.7d0)
    call expand_array(amax = 1.d0/(1.d0 + max(this%zrecomb, 1.d0)), stepfac = 0.96d0, stepmin = step_ini*0.5d0)
    call expand_array(amax = 1.d0/(1.d0 + max(this%zrecomb-200.d0, 1.d0)), stepfac = 1.02d0, stepmax = step_ini*0.8d0)    
    call expand_array(amax = 1.d0/(1.d0 + max(this%zrecomb_end-100.d0, 1.d0)), stepfac = 1.02d0, stepmax = step_ini*1.5d0)    
    if(this%optre.gt.0.01d0)then
       call expand_array(amax =  min(1.d0/(1.d0 + this%zre + this%deltaz*15.d0), 0.5d0), stepfac = 1.03d0, stepmax = step_ini*10.d0)
       call expand_array(amax =  min(1.d0/(1.d0 + this%zre), 0.5d0), stepfac = 0.97d0, stepmin = step_ini*2.d0)
       call expand_array(amax =  min(1.d0/(1.d0 + this%zre-this%deltaz*4.d0), 0.5d0), stepfac = 1.03d0, stepmax = step_ini*5.d0)
       call expand_array(amax = 0.75d0, stepfac = 1.03d0, stepmax = step_ini*15.d0)
    else
       call expand_array(amax = 0.75d0, stepfac = 1.02d0, stepmax = step_ini*15.d0)
    endif
    call expand_array(amax = 0.99999d0, stepfac = 0.96d0, stepmin = step_ini*2.d0)
    source%ntau = n
    if(allocated(source%a))then
       if(size(source%a).ne.n)then
          deallocate(source%a, source%tau, source%chi, source%dtau, source%tauc, source%lna, source%omega_rad, source%vis)
          allocate(source%a(n), source%tau(n), source%dtau(n), source%chi(n), source%tauc(n), source%lna(n), source%omega_rad(n), source%vis(n), source%omega_de(n))
       endif
    else
       allocate(source%a(n),source%tau(n), source%dtau(n), source%chi(n),  source%tauc(n), source%lna(n), source%omega_rad(n), source%vis(n), source%omega_de(n))
    endif
    source%a = a(1:n)
    source%tau = tau(1:n)
    source%dtau(1) = source%tau(2) - source%tau(1)
    do i = 2, n-1
       source%dtau(i) = (source%tau(i+1) - source%tau(i-1))/2.d0
    enddo
    source%dtau(n) = source%tau(n) - source%tau(n-1)
    source%chi = this%tau0 - source%tau
    source%lna = log(source%a)
    !$omp parallel do
    do i=1, n
       source%tauc(i) = this%taucofa(source%a(i))
       source%vis(i) = this%visofa(source%a(i))
       source%omega_rad(i) = this%Omega_r/( this%Omega_r + source%a(i)*(this%Omega_m+(1.d0-this%Omega_m)*source%a(i)**3))  !!this is just an approximation, only used for optimization (drop radiation when omega_rad is very small)
       source%omega_de(i) = O0_DE(this)%rhoa4(source%a(i))/this%rhoa4(source%a(i))
    enddo
    !$omp end parallel do
    source%index_vis_max = coop_maxloc(source%vis)
    source%index_vis_end = source%index_vis_max
    viscut = source%vis(source%index_vis_max)*0.003d0
    do while(source%index_vis_end .lt. n .and. source%vis(source%index_vis_end) .gt. viscut )
       source%index_vis_end = source%index_vis_end + 1
    enddo
    viscut = source%vis(source%index_vis_max)*1.d-4
    source%index_vis_start = source%index_vis_max    
    do while(source%index_vis_start .gt. 1 .and. source%vis(source%index_vis_start) .gt. viscut )
       source%index_vis_start = source%index_vis_start - 1
    enddo
  contains
    subroutine expand_array(amax, stepfac, stepmin, stepmax)
      COOP_REAL::amax
      COOP_REAL,optional::stepfac, stepmin, stepmax
      do
         n =n+1
         if(n.gt.nmax) call coop_return_error("set_source_tau", "nmax<n", "stop")
         tau(n) = tau(n-1) + step*step_factor
         a(n) = this%aoftau(tau(n))
         if(a(n).gt. amax)then
            n = n-1
            return
         endif
         if(present(stepfac))then
            step = step * stepfac
            if(present(stepmin))then
               step = max(step, stepmin)
            endif
            if(present(stepmax))then
               step = min(step, stepmax)
            endif
         endif
      enddo
    end subroutine expand_array
    
  end subroutine coop_cosmology_firstorder_set_source_tau


  subroutine coop_cosmology_firstorder_set_source_k(this, source, n, weight)
    class(coop_cosmology_firstorder)::this
    class(coop_cosmology_firstorder_source)::source
    COOP_INT n, i, iq, j
    COOP_REAL weight, dlnk, tauc_cut, omega_rad_cut
    this%k_pivot = this%kMpc_pivot/this%H0Mpc()
    source%kweight = weight
    source%nk = n
    if(allocated(source%k))then
       if(size(source%k).ne.n)then
          deallocate(source%k, source%dk, source%index_tc_off, source%kop)
          allocate(source%k(n), source%dk(n), source%index_tc_off(n),  source%kop(n))
       endif
       if(size(source%k_dense, 2).ne.n .or. size(source%k_dense, 1).ne. coop_k_dense_fac )then
          deallocate(source%k_dense, source%ws_dense, source%wt_dense, source%ps_dense)
          allocate(source%k_dense(coop_k_dense_fac, n), source%ws_dense(coop_k_dense_fac, n), source%wt_dense(coop_k_dense_fac, n), source%ps_dense(coop_k_dense_fac, n) )
       endif
    else
       allocate(source%k(n), source%dk(n), source%index_tc_off(n),  source%kop(n))
       allocate(source%k_dense(coop_k_dense_fac, n), source%ws_dense(coop_k_dense_fac, n), source%wt_dense(coop_k_dense_fac, n), source%ps_dense(coop_k_dense_fac, n) )
    endif
    do i=1, coop_k_dense_fac
       source%a_dense(i) = dble(i)/coop_k_dense_fac
       source%b_dense(i) =  1.d0 - source%a_dense(i)
       source%a2_dense(i) = source%a_dense(i)*(source%a_dense(i)**2-1.d0)
       source%b2_dense(i) = source%b_dense(i)*(source%b_dense(i)**2-1.d0)
    enddo

    if(coop_force_kminkmax)then
       source%kmin = coop_power_kmin
       source%kmax = coop_power_kmax
    else
       source%kmin = 0.4d0/this%distlss
#if DO_ZETA_TRANS
       source%kmax = (min(max(1500, coop_Cls_lmax(source%m)), 3500)*2.5d0)/this%distlss
#else	 
       source%kmax = (min(max(1500, coop_Cls_lmax(source%m)), 3500)*2.2d0)/this%distlss
#endif
    endif
    call source%k2kop(source%kmin, source%kopmin)
    call source%k2kop(source%kmax, source%kopmax)
    source%dkop = (source%kopmax-source%kopmin)/(n-1)
    source%dkop_dense = source%dkop/coop_k_dense_fac    
    call coop_set_uniform(n, source%kop, source%kopmin, source%kopmax)

    !$omp parallel do private(i, j, dlnk)
    do i=1, n
       call source%kop2k(source%kop(i), source%k(i), source%dkop,  source%dk(i))
       do j=1, coop_k_dense_fac
          call source%kop2k(source%kop(i)+(j-coop_k_dense_fac)*source%dkop_dense, source%k_dense(j,i), source%dkop_dense,  dlnk)
          dlnk = dlnk/source%k_dense(j, i)
!          source%ps_dense(j, i) = this%psofk(source%k_dense(j, i))
          source%ps_dense(j, i) = this%pssmooth%eval(source%k_dense(j, i))
          source%ws_dense(j, i) = dlnk*source%ps_dense(j, i)
!          source%wt_dense(j, i) = this%ptofk(source%k_dense(j, i))*dlnk
          source%wt_dense(j, i) = this%ptsmooth%eval(source%k_dense(j, i))*dlnk
       enddo
    enddo
    !$omp end parallel do

    !!set index_tc_off, index_de_perturb_on
    if(source%ntau .gt. 0 .and. allocated(source%tauc))then
       !$omp parallel do private(tauc_cut)
       do i = 1, n
          tauc_cut = min(0.05, coop_cosmology_firstorder_tc_cutoff*(1.d0+source%k(i)/500.d0))/source%k(i)  !!for very large k we use less accurate scheme to speed up the code
          source%index_tc_off(i) = max(1, min(coop_left_index(source%index_vis_max, source%tauc(1:source%index_vis_max), tauc_cut), coop_left_index(source%index_vis_max, source%tauc(1:source%index_vis_max)/source%tau(1:source%index_vis_max), coop_cosmology_firstorder_tc_cutoff), source%ntau))
       enddo
       !$omp end parallel do
       !!allow early DE?
       !!source%index_de_perturb_on = coop_right_index(source%ntau, source%omega_de, 1.d-6)
       source%index_de_perturb_on = max(maxval(source%index_tc_off)+1, coop_right_index(source%ntau, source%omega_de, 1.d-6) )
    else
       call coop_return_error("set_source_k", "You need to call set_source_tau before calling set_source_k", "stop")
    endif
    
    !!set index_massivenu_on
    if(this%mnu_by_Tnu*source%a(source%ntau) .lt. coop_pert_massivenu_threshold(1))then
       source%index_massivenu_on = source%ntau + 1
       source%index_massivenu_cold = source%ntau + 1
       return
    endif
    source%index_massivenu_cold = source%ntau + 1
    do while(this%mnu_by_Tnu*source%a(source%index_massivenu_cold-1) .gt. coop_pert_massivenu_cold_threshold)
       source%index_massivenu_cold = source%index_massivenu_cold - 1
       if(source%index_massivenu_cold .le. 1) exit
    enddo

    source%index_massivenu_on(coop_pert_default_nq) = source%index_massivenu_cold
    iq = coop_pert_default_nq
    if(source%index_massivenu_on(iq) .gt. 1)then
       do while(this%mnu_by_Tnu*source%a(source%index_massivenu_on(iq)-1) .gt. coop_pert_massivenu_threshold(iq))
          source%index_massivenu_on(iq) = source%index_massivenu_on(iq)-1
          if(source%index_massivenu_on(iq) .le. 1)exit
       enddo
    endif
    do iq = coop_pert_default_nq - 1, 1, -1
       source%index_massivenu_on(iq) =  source%index_massivenu_on(iq+1)
       if(source%index_massivenu_on(iq) .le. 1)cycle
       do while(this%mnu_by_Tnu*source%a(source%index_massivenu_on(iq)-1) .gt. coop_pert_massivenu_threshold(iq))
          source%index_massivenu_on(iq) = source%index_massivenu_on(iq)-1
          if(source%index_massivenu_on(iq) .le. 1)exit
       enddo
    enddo
  end subroutine coop_cosmology_firstorder_set_source_k



   function coop_cosmology_firstorder_source_interpolate(source, k, chi) result(s)
    COOP_INT ik, ichi
    COOP_INT ilow, iup
    class(coop_cosmology_firstorder_source)::source
    COOP_REAL k, chi, kop, rk, a, b, rchi
    COOP_REAL s(source%nsrc)
    if(chi .le. source%chi(source%ntau) .or. chi .ge. source%chi(1))then
       s = 0.d0
       return
    endif
    call source%k2kop(k, kop)
    rk = (kop - source%kopmin)/source%dkop + 1.d0
    ik = floor(rk)
    if(ik.le. 0 .or. ik.ge.source%nk)then
       s = 0.d0
       return
    endif
    a = rk - ik
    b = 1.d0 - a
    ilow = 1
    iup = source%ntau
    do while(iup - ilow .gt. 1)
       ichi = (iup + ilow)/2
       if(source%chi(ichi) .gt. chi)then
          ilow = ichi
       else
          iup = ichi
       endif
    enddo
    rchi = source%chi(ilow) - chi
    s = ((source%s(:, ik, ilow)+source%s2(:, ik, ilow)*(b**2-1.d0))*b + (source%s(:, ik+1, ilow) + source%s2(:, ik+1, ilow)*(a**2-1.d0))*a)*(1.d0-rchi) &
         +  ((source%s(:, ik, iup)+source%s2(:, ik, iup)*(b**2-1.d0))*b + (source%s(:, ik+1, iup) + source%s2(:, ik+1, iup)*(a**2-1.d0))*a)*rchi
  end function coop_cosmology_firstorder_source_interpolate


  function coop_cosmology_firstorder_source_interpolate_one(source, k, chi, isrc) result(s)
    COOP_INT ik, ichi, isrc
    COOP_INT ilow, iup
    class(coop_cosmology_firstorder_source)::source
    COOP_REAL k, chi, kop, rk, a, b, rchi
    COOP_REAL s
    if(chi .le. source%chi(source%ntau) .or. chi .ge. source%chi(1))then
       s = 0.d0
       return
    endif
    call source%k2kop(k, kop)
    rk = (kop - source%kopmin)/source%dkop + 1.d0
    ik = floor(rk)
    if(ik.le. 0 .or. ik.ge.source%nk)then
       s = 0.d0
       return
    endif
    a = rk - ik
    b = 1.d0 - a
    ilow = 1
    iup = source%ntau
    do while(iup - ilow .gt. 1)
       ichi = (iup + ilow)/2
       if(source%chi(ichi) .gt. chi)then
          ilow = ichi
       else
          iup = ichi
       endif
    enddo
    rchi = source%chi(ilow) - chi
    s = ((source%s(isrc, ik, ilow)+source%s2(isrc, ik, ilow)*(b**2-1.d0))*b + (source%s(isrc, ik+1, ilow) + source%s2(isrc, ik+1, ilow)*(a**2-1.d0))*a)*(1.d0-rchi) &
         +  ((source%s(isrc, ik, iup)+source%s2(isrc, ik, iup)*(b**2-1.d0))*b + (source%s(isrc, ik+1, iup) + source%s2(isrc, ik+1, iup)*(a**2-1.d0))*a)*rchi
  end function coop_cosmology_firstorder_source_interpolate_one


  
  subroutine coop_cosmology_firstorder_source_get_transfer(source, l, trans)
    class(coop_cosmology_firstorder_source)::source
    COOP_INT,parameter::lc1 = 250, lc2 = 280
    COOP_INT::l, limber_start
    COOP_REAL,dimension(:,:,:)::trans
    COOP_INT::n, ik, idense, itau, ikmin, ikmax, limber_ikmin, itau_cut
    COOP_REAL::jl, xmin, chi_source, wl, wr, const,  kchicut, x!, dampfac

    if(size(trans,2).ne. coop_k_dense_fac .or. size(trans, 3).ne. source%nk .or. size(trans, 1) .ne. source%nsrc) call coop_return_error("get_transfer", "wrong size", "stop")
    select case(source%m)
    case(0)
       limber_start = coop_index_source_Len
    case(2)
       limber_start = coop_index_source_B + 1       
    end select
    trans = 0.d0
    call coop_jl_check_init(l)
    call coop_jl_startpoint(l, xmin)
    ikmin = max( coop_right_index(source%nk, source%k, xmin/source%chi(source%index_vis_start)),  2 )
    ikmax = min( source%nk,  ceiling(max(l*5.d0, 800.d0)/source%distlss) )
  !  dampfac = 1.2d0/source%k(ikmax)
    !!l < 20;  Want reionization; brute-force integral
    if(l .lt. lc2)then 
       !$omp parallel do private(ik, itau, idense, jl, x)
       do ik=ikmin, ikmax
          do idense = 1, coop_k_dense_fac 
             do itau = source%index_vis_start, source%ntau
                x = source%k_dense(idense, ik)*source%chi(itau)
                if(x .lt. xmin)exit
                jl = coop_jl(l, x )
                trans(:, idense, ik) = trans(:, idense, ik) + (jl* source%dtau(itau)) * COOP_INTERP_SOURCE(source, :, idense, ik, itau)
             enddo
          enddo
       enddo
       !$omp end parallel do
       if(source%m.eq.0 .and. l .gt. lc1)then
          trans(limber_start:source%nsrc, :,:) = trans(limber_start:source%nsrc, :,:) * cos(coop_pio2*(l-lc1)/(lc2-lc1))**2
          const = sqrt(coop_pi/2.d0/l)* sin(coop_pio2*(l-lc1)/(lc2-lc1))**2
          !$omp parallel do private(ik, itau, idense, chi_source, wl, wr)       
          do ik = ikmin, ikmax
             do idense = 1, coop_k_dense_fac
                !limber contribution
                chi_source = l/source%k_dense(idense, ik)
                itau = coop_left_index(source%ntau, source%chi, chi_source)
                if(itau .ge. 1 .and. itau .lt. source%ntau)then
                   wr = (source%chi(itau) - chi_source)/(source%chi(itau) - source%chi(itau + 1))
                   wl = 1.d0 - wr
                   trans(limber_start:source%nsrc, idense, ik) =  trans(limber_start:source%nsrc, idense, ik) + const/source%k_dense(idense, ik) *  (( COOP_INTERP_SOURCE(source, limber_start:source%nsrc, idense, ik, itau) *  wl + COOP_INTERP_SOURCE(source, limber_start:source%nsrc, idense, ik, itau+1) *  wr ))
                endif
             enddo
          enddo
          !$omp end parallel do
       endif
       
       goto 100
    endif
    const = sqrt(coop_pi/2.d0/l)        
    if(l .lt. coop_limber_ell)then
       !$omp parallel do private(ik, itau, idense, jl, x)
       do ik=ikmin, ikmax
          do idense = 1, coop_k_dense_fac 
             do itau = source%index_vis_start, source%ntau
                x = source%k_dense(idense, ik)*source%chi(itau)
                if(x .lt. xmin)exit
                jl = coop_jl(l, x )
                trans(1:limber_start-1, idense, ik) = trans(1:limber_start-1, idense, ik) + (jl* source%dtau(itau)) * COOP_INTERP_SOURCE(source, 1:limber_start-1, idense, ik, itau) 
             enddo
          enddo
       enddo
       !$omp end parallel do
       if(source%m.eq.0)then
          !$omp parallel do private(ik, itau, idense, chi_source, wl, wr)       
          do ik = ikmin, ikmax
             do idense = 1, coop_k_dense_fac
                !limber contribution
                chi_source = l/source%k_dense(idense, ik)
                itau = coop_left_index(source%ntau, source%chi, chi_source)
                if(itau .ge. 1 .and. itau .lt. source%ntau)then
                   wr = (source%chi(itau) - chi_source)/(source%chi(itau) - source%chi(itau + 1))
                   wl = 1.d0 - wr
                   trans(limber_start:source%nsrc, idense, ik) =  trans(limber_start:source%nsrc, idense, ik) + const/source%k_dense(idense, ik) *  (( COOP_INTERP_SOURCE(source, limber_start:source%nsrc, idense, ik, itau) *  wl + COOP_INTERP_SOURCE(source, limber_start:source%nsrc, idense, ik, itau+1) *  wr ))
                endif
             enddo
          enddo
          !$omp end parallel do
       endif
       goto 100
    endif
    !!otherwise it is a mixture of limber contribution and the contribution from last scattering surface:
    !$omp parallel do private(ik, itau, idense, chi_source, wl, wr, x, jl)
    do ik = ikmin, ikmax
       do idense = 1, coop_k_dense_fac
          !limber contribution
          chi_source = l/source%k_dense(idense, ik)
          itau = coop_left_index(source%ntau, source%chi, chi_source)
          if(itau .ge. 1 .and. itau .lt. source%ntau)then
             wr = (source%chi(itau) - chi_source)/(source%chi(itau) - source%chi(itau + 1))
             wl = 1.d0 - wr
             if(itau .gt. source%index_vis_end)then
                trans(:, idense, ik) =  trans(:, idense, ik) + const/source%k_dense(idense, ik) *  ( COOP_INTERP_SOURCE(source, :, idense, ik, itau) *  wl + COOP_INTERP_SOURCE(source, :, idense, ik, itau+1) *  wr )
             elseif(source%m.eq.0)then
                trans(limber_start:source%nsrc, idense, ik) =  trans(limber_start:source%nsrc, idense, ik) + const/source%k_dense(idense, ik) *  (( COOP_INTERP_SOURCE(source, limber_start:source%nsrc, idense, ik, itau) *  wl + COOP_INTERP_SOURCE(source, limber_start:source%nsrc, idense, ik, itau+1) *  wr ))            
             endif
          endif
          !!recombination contribution
          do itau = source%index_vis_start, source%index_vis_end
             x = source%k_dense(idense, ik)*source%chi(itau)
             if(x .lt. xmin)exit
             jl = coop_jl(l, x )
             trans(1:limber_start-1, idense, ik) = trans(1:limber_start-1, idense, ik) + (jl* source%dtau(itau)) * COOP_INTERP_SOURCE(source, 1:limber_start-1, idense, ik, itau) 

          enddo
       enddo
    enddo
    !$omp end parallel do
100 continue
#if DO_ZETA_TRANS
    if(coop_zeta_single_slice_chi .gt. 0.d0 .and. .not. coop_zeta_user_specified_weight%initialized)then
       !$omp parallel do private(ik , idense)
       do ik = 1, source%nk
          do idense = 1, coop_k_dense_fac
             trans(coop_index_source_zeta, idense, ik) = -coop_jl(l, source%k_dense(idense, ik)*coop_zeta_single_slice_chi)
          enddo
       enddo
       !$omp end parallel do
    endif
#endif    
  end subroutine coop_cosmology_firstorder_source_get_transfer

  subroutine coop_cosmology_firstorder_source_get_Cls_limber(source, l, Cls)
    class(coop_cosmology_firstorder_source)::source
    COOP_INT l, ik, idense, itau
    COOP_REAL,dimension(coop_num_cls),intent(OUT)::Cls
    COOP_REAL::src(source%nsrc), kop, rk
    Cls = 0.d0
    do itau = source%index_vis_end + 1, source%ntau
       call source%k2kop(l/source%chi(itau), kop)
       rk = (kop-source%kopmin)/source%dkop+1.d0
       ik = ceiling(rk)
       if(ik.lt.2 .or. ik .gt. source%nk)cycle
       idense = coop_k_dense_fac - nint((ik - rk)*coop_k_dense_fac)
       if(idense .eq. 0)then
          ik = ik - 1
          if(ik.eq.1)cycle
          idense = coop_k_dense_fac
       endif
       src = COOP_INTERP_SOURCE(source, :, idense, ik, itau)
       rk = source%ps_dense(idense,ik)*source%chi(itau)*source%dtau(itau)
       Cls(coop_index_ClTT) = Cls(coop_index_ClTT) + rk*src(coop_index_source_T)**2
       Cls(coop_index_ClEE) = Cls(coop_index_ClEE) + rk*src(coop_index_source_E)**2
       Cls(coop_index_ClTE) = Cls(coop_index_ClTE) + rk*src(coop_index_source_T)*src(coop_index_source_E)
       if(source%m .eq. 0)then
          if(source%nsrc .ge. 3)then
             Cls(coop_index_ClLenLen) = Cls(coop_index_ClLenLen) + rk*src(coop_index_source_Len)**2
             Cls(coop_index_ClTLen) = Cls(coop_index_ClTLen) + rk*src(coop_index_source_T)*src(coop_index_source_Len)
#if DO_ZETA_TRANS             
             if(source%nsrc .ge. 4)then
                Cls(coop_index_ClTzeta) = Cls(coop_index_ClTzeta) + rk*src(coop_index_source_zeta)*src(coop_index_source_T)
                Cls(coop_index_ClEzeta) = Cls(coop_index_ClTzeta) + rk*src(coop_index_source_zeta)*src(coop_index_source_E)
                Cls(coop_index_Clzetazeta) = Cls(coop_index_ClTzeta) + rk*src(coop_index_source_zeta)**2
             endif
#endif             
          endif
       else
          Cls(coop_index_ClBB) = Cls(coop_index_ClBB) + rk*src(coop_index_source_B)**2
       endif
    enddo
    if(source%m .eq. 0 )then
       do itau = 1, source%index_vis_end
          call source%k2kop(l/source%chi(itau), kop)
          rk = (kop-source%kopmin)/source%dkop+1.d0
          ik = ceiling(rk)
          if(ik.lt.2 .or. ik .gt. source%nk)cycle
          idense = coop_k_dense_fac - nint((ik - rk)*coop_k_dense_fac)
          if(idense .eq. 0)then
             ik = ik - 1
             if(ik.eq.1)cycle
             idense = coop_k_dense_fac
          endif
          src = COOP_INTERP_SOURCE(source, :, idense, ik, itau)
          rk = source%ps_dense(idense,ik)*source%chi(itau)*source%dtau(itau)

          if(source%nsrc .ge. 3)then
             Cls(coop_index_ClLenLen) = Cls(coop_index_ClLenLen) + rk*src(coop_index_source_Len)**2
             Cls(coop_index_ClTLen) = Cls(coop_index_ClTLen) + rk*src(coop_index_source_T)*src(coop_index_source_Len)
#if DO_ZETA_TRANS             
             if(source%nsrc .ge. 4)then
                Cls(coop_index_ClTzeta) = Cls(coop_index_ClTzeta) + rk*src(coop_index_source_zeta)*src(coop_index_source_T)
                Cls(coop_index_ClEzeta) = Cls(coop_index_ClTzeta) + rk*src(coop_index_source_zeta)*src(coop_index_source_E)
                Cls(coop_index_Clzetazeta) = Cls(coop_index_ClTzeta) + rk*src(coop_index_source_zeta)**2
             endif
#endif             
          endif
       enddo
    endif
    Cls = Cls*((2.d0*coop_pi**2)/dble(l)**3)
  end subroutine coop_cosmology_firstorder_source_get_Cls_limber


