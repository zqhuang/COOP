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
    !!this computes 4 \pi \int_0^\infty |j_l(kr)|^2 (k^3P(k)/(2\pi^2)) d\ln k
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
    COOP_REAL,parameter::k_trunc = 10.d0  !!for tiny k the integrated-by-part term is negligible.
    class(coop_cosmology_firstorder_source)::source
    COOP_REAL s2(source%ntau), sum1
    COOP_INT ik, i
    if(source%k(ik).lt. k_trunc/5.d0)return
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
    if(allocated(this%k))deallocate(this%k, this%dk, this%index_tc_off, this%kop,this%index_rad_off)
    this%nk = 0
    if(allocated(this%tau))deallocate(this%tau, this%chi, this%dtau, this%a, this%tauc, this%lna, this%omega_rad, this%vis)
    this%ntau = 0
    if(allocated(this%s))deallocate(this%s, this%s2, this%saux)
    if(allocated(this%k_dense))deallocate(this%k_dense, this%ws_dense, this%wt_dense, this%ps_dense)
  end subroutine coop_cosmology_firstorder_source_free

  subroutine coop_cosmology_firstorder_free(this)
    class(coop_cosmology_firstorder)::this
    COOP_INT m, i
    do m=0,2
       call this%source(m)%free()
    enddo
    call this%Ps%free()
    call this%Pt%free()
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
  end subroutine coop_cosmology_firstorder_free


  subroutine coop_cosmology_firstorder_set_power(this, power, args)
!!the format of primordial power
!!subroutine power(k/k_pivot, ps, pt, cosmology, args)
!!input kMpc and args, output ps and pt    
    integer,parameter::n = 4096
    class(coop_cosmology_firstorder)::this
    external power
    type(coop_arguments) args
    COOP_REAL k(n), ps(n), pt(n)
    COOP_INT i
    this%k_pivot = this%kMpc_pivot/this%H0Mpc()
    call coop_set_uniform(n, k, coop_power_lnk_min, coop_power_lnk_max)
    k = exp(k)
    do i=1, n
       call power(k(i)/this%k_pivot, ps(i), pt(i), this, args)
    end do
    call this%ps%init(n = n, xmin = k(1), xmax = k(n), f = ps, xlog = .true., ylog = .true., fleft = ps(1), fright = ps(n), check_boundary = .false.)
    if(any(pt.eq.0.d0))then
       call this%pt%init(n = n, xmin = k(1), xmax = k(n), f = pt, xlog = .true., ylog = .false., fleft = pt(1), fright = pt(n), check_boundary = .false.)
    else
       call this%pt%init(n = n, xmin = k(1), xmax = k(n), f = pt, xlog = .true., ylog = .true., fleft = pt(1), fright = pt(n), check_boundary = .false.)
    endif
    this%As = this%psofk(this%k_pivot)
    this%ns = this%ps%derivative_bare(log(this%k_pivot)) + 1.d0
    this%nrun = this%ps%derivative2_bare(log(this%k_pivot))
    this%r =  this%ptofk(this%k_pivot)/this%As
    this%has_tensor = any(pt .gt. 1.d-15)
    if(this%r .eq. 0.d0)then
       this%nt = 0.d0
    else
       this%nt = this%pt%derivative(this%k_pivot)*this%k_pivot/(this%r*this%As)
    endif
  end subroutine coop_cosmology_firstorder_set_power

  subroutine coop_cosmology_firstorder_set_xe(this)
    class(coop_cosmology_firstorder)::this
    integer i, m
    integer, parameter::n = 8192
    COOP_REAL ekappa(n), a(n), dkappadinva(n), dkappadtau(n), vis(n)
    call this%setup_background()
    call this%set_klms()
    this%index_baryon = this%index_of("Baryon", .true.)
    this%index_radiation = this%index_of("Radiation", .true.)
    this%index_cdm = this%index_of("CDM", .true.)
    this%index_Nu = this%index_of("Massless Neutrinos", .true.)
    this%index_massiveNu = this%index_of("Massive Neutrinos")
    this%index_de = this%index_of("Dark Energy", .true.)
    this%Omega_b = O0_BARYON(this)%Omega
    this%ombh2 = this%Omega_b * this%h()**2    
    this%Omega_c = O0_CDM(this)%Omega
    this%omch2 = this%Omega_c*O0_CDM(this)%rhoa3_ratio(1.d-10)*this%h()**2
    this%Omega_g = O0_RADIATION(this)%Omega 
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

    this%Omega_m = this%Omega_b + this%Omega_c + this%Omega_massivenu 
    this%Omega_r =  this%Omega_g + this%Omega_nu 
    this%bbks_keq =  (0.073*coop_SI_c/1.d5) * this%h_value * this%omega_m * exp(- this%omega_b - sqrt(2.*this%h_value)*this%omega_b/this%omega_m)

    this%tau_eq = this%conformal_time(this%a_eq)

    this%dkappadtau_coef = this%species(this%index_baryon)%Omega * this%h() * coop_SI_sigma_thomson * (coop_SI_rhocritbyh2/coop_SI_c**2) * coop_SI_hbyH0 * coop_SI_c/ coop_SI_m_H * (1.d0 - this%YHe())
    if(this%do_reionization)then
       this%reionFrac = 1.d0 + this%YHe()/(coop_m_He_by_m_H * (1.d0-  this%YHe()))
    else
       this%reionFrac = 0.d0
    endif

    call this%set_zre_from_optre()
    call coop_recfast_get_xe(this, this%xe, this%Tb, this%reionFrac, this%zre, this%deltaz)
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
    do while(this%vis%eval(1.d0/(this%zrecomb_start+1.d0))/this%maxvis .gt. 1.d-4 .and. this%zrecomb_start .lt. 2.d4)
       this%zrecomb_start = this%zrecomb_start + 30.d0
    enddo

    this%zrecomb_end = this%zrecomb - 10.d0
    do while(this%vis%eval(1.d0/(this%zrecomb_end+1.d0))/this%maxvis .gt. 1.d-3 .and. this%zrecomb_end .gt. 100.d0)
       this%zrecomb_end = this%zrecomb_end  - 10.d0      
    enddo

    this%distlss = this%comoving_distance(1.d0/(1.d0+this%zrecomb))
    this%tau0 = this%conformal_time(coop_scale_factor_today)
    
  end subroutine coop_cosmology_firstorder_set_xe

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

  subroutine coop_cosmology_firstorder_standard_power(kbykpiv, ps, pt, cosmology, args)
    type(coop_cosmology_firstorder)::cosmology
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
    COOP_REAL zstar, theta, astar, rs, da, ombh2
    zstar = coop_zrecomb_fitting(this%ombh2, this%omch2)
    astar = 1.d0/(1.d0+zstar)
    rs = coop_integrate(dsoundda, 1.d-8, astar, 1.d-7)  !!to be consistent with CosmoMC
    da = coop_r_of_chi(coop_integrate(dtauda, astar, 1.d0, 1.d-7), this%Omega_k())
    theta  = rs/da
  contains

    function dsoundda(a) result(dsda)  !!use approximations, to be consistent with CosmoMC
      COOP_REAL a, dsda, R
      R = this%ombh2 * 3.d4*a
      dsda = dtauda(a)/sqrt(3.d0*(1.d0+R))
!      dsda = (1.d0/coop_sqrt3) / sqrt(1.d0 + 0.75d0/ this%Omega_radiation()/this%h()**2 * a * ombh2 )/ this%Hasq(a)
    end function dsoundda

    function dtauda(a) 
      COOP_REAL a, dtauda
      dtauda = 1.d0/this%Hasq(a)
    end function dtauda
  end function coop_cosmology_firstorder_cosmomc_theta
	    
