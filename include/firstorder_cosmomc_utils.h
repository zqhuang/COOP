  subroutine coop_cosmology_firstorder_allocate_source(this, m, source, tau, k, indices)
    class(coop_cosmology_firstorder)::this
    class(coop_cosmology_firstorder_source)::source
    COOP_INT :: m
    COOP_REAL,dimension(:),optional::tau
    COOP_REAL,dimension(:), optional::k
    COOP_INT,dimension(:), optional::indices
    source%distlss = this%distlss
    source%m = m
    source%nsrc = coop_num_sources(m)
    source%nsaux = coop_num_saux(m)
    if(present(tau))then
       if(present(indices))then
          if(size(indices).eq.size(tau))then
             call this%set_source_tau(source, coop_source_tau_step_factor(m), tau_wanted = tau, indices = indices)
          else
             stop "allocate_source: the size of tau and indices must be the same"
          endif
       else
          call this%set_source_tau(source, coop_source_tau_step_factor(m), tau_wanted = tau)
       endif
    else
       call this%set_source_tau(source, coop_source_tau_step_factor(m))
    endif
    if(present(k))then
       call this%set_source_given_k(source, k)
    else
       call this%set_source_k(source, coop_source_k_n(m), coop_source_k_weight(m))
    endif
    if(allocated(source%s))then
       if(size(source%s, 1) .ne. source%ntau .or. size(source%s, 2) .ne. source%nk)then
          deallocate(source%s, source%s2, source%saux)
          allocate(source%s(source%nsrc, source%nk , source%ntau), source%s2(source%nsrc, source%nk, source%ntau), source%saux(source%nsaux, source%nk, source%ntau) )
       endif
    else
       allocate(source%s(source%nsrc, source%nk , source%ntau), source%s2(source%nsrc, source%nk, source%ntau), source%saux(source%nsaux, source%nk, source%ntau) )
    endif
  end subroutine coop_cosmology_firstorder_allocate_source
  

  subroutine coop_cosmology_firstorder_set_source_given_tau(this, source, tau)
    class(coop_cosmology_firstorder_source)::source
    class(coop_cosmology_firstorder)::this
    COOP_REAL,dimension(:)::tau
    COOP_INT i, n
    n=size(tau)
    source%ntau = n
    if(allocated(source%a))then
       if(size(source%a).ne.n)then
          deallocate(source%a, source%tau, source%chi, source%dtau, source%tauc, source%lna)
          allocate(source%a(n), source%tau(n), source%dtau(n), source%chi(n), source%tauc(n), source%lna(n))
       endif
    else
       allocate(source%a(n), source%tau(n), source%dtau(n), source%chi(n),  source%tauc(n), source%lna(n))
    endif

    source%tau = tau
    !$omp parallel do
    do i=1, n
       source%a(i) = this%aoftau(source%tau(i))
    enddo
    !$omp end parallel do

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
    enddo
    !$omp end parallel do
  end subroutine coop_cosmology_firstorder_set_source_given_tau

  subroutine coop_cosmology_firstorder_set_source_given_k(this, source, k)
    class(coop_cosmology_firstorder)::this
    class(coop_cosmology_firstorder_source)::source
    COOP_REAL,dimension(:):: k
    COOP_INT :: n, i, iq
    this%k_pivot = this%kMpc_pivot/this%H0Mpc()
    source%kweight = 1.d0
    n = size(k)
    source%nk = n
    if(allocated(source%k))then
       if(size(source%k).ne.n)then
          deallocate(source%k, source%dk, source%index_tc_off, source%kop)
          allocate(source%k(n), source%dk(n), source%index_tc_off(n), source%kop(n))
       endif
    else
       allocate(source%k(n), source%dk(n), source%index_tc_off(n), source%kop(n))
    endif
    source%kmin = k(1)
    source%kmax = k(n)
    source%k = k

    !!set index_tc_off
    if(source%ntau .gt. 0 .and. allocated(source%tauc))then
       source%index_tc_off = 1
       do i = 1, n
          do while(source%index_tc_off(i) .lt. source%ntau - 1 .and. source%tauc(source%index_tc_off(i)+1)*source%k(i) .le. coop_cosmology_firstorder_tc_cutoff)
             source%index_tc_off(i)= source%index_tc_off(i)+1
          enddo
       enddo
    else
       call coop_return_error("set_source_given_k", "You need to call set_source_tau before calling set_source_given_k", "stop")
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
  end subroutine coop_cosmology_firstorder_set_source_given_k
						    
  subroutine coop_cosmology_firstorder_camb_DoSourceK(this, m, ik, kMpc, tauMpc, source, num_trans, tauMpc_trans, trans)
    class(coop_cosmology_firstorder)::this
    type(coop_cosmology_firstorder_source)::s
    COOP_INT m, ik
    COOP_INT, optional::num_trans
    COOP_REAL,dimension(:),optional::tauMpc_trans
    COOP_REAL kMpc, tauMpc(:), h0mpc, psi, phinewt, phiweyl, vc, h0tau
    COOP_REAL,dimension(:,:,:)::source
    COOP_SINGLE,dimension(:,:,:),optional::trans
    COOP_INT::ntau, i, itf
    COOP_INT,dimension(:),allocatable::indices
    COOP_REAL::  a, kbyH0
    ntau = size(tauMpc)
    allocate(indices(ntau))
    h0mpc = this%H0Mpc()
    kbyH0 = kMpc/h0mpc
    call this%allocate_source(m = m, source = s, k = (/ kbyH0 /), tau = tauMpc*h0mpc, indices=indices)
    call this%compute_source_k(s, 1)
    source(ik, 1, :) = s%s(1, 1, indices)*h0mpc
    source(ik, 2, :) = s%s(2, 1, indices)*h0mpc 
    source(ik, 3, :) = s%s(3, 1, indices)*h0mpc
    if(m .eq. 0 .and. present(trans))then
       trans(1, ik, :) = kMpc/this%h()       
       do itf = 1, num_trans
          h0tau = tauMpc_trans(itf)*h0mpc
          call coop_linear_interp(s%ntau, s%tau, s%saux(3, 1, :), h0tau, psi)
          call coop_linear_interp(s%ntau, s%tau, s%saux(2, 1, :), h0tau, phiweyl)
          call coop_linear_interp(s%ntau, s%tau, s%saux(6, 1, :), h0tau, vc)          
          a = this%aoftau(h0tau)
          phinewt = phiweyl - psi
          if(this%index_massivenu.ne.0)then
             trans(7, ik, itf) = phinewt/(1.5d0*h0mpc**2*(this%Omega_b*O0_BARYON(this)%density_ratio(a)+this%Omega_c*O0_CDM(this)%density_ratio(a) + this%Omega_massivenu*O0_MASSIVENU(this)%density_ratio(a) )*a**2)
          else
             trans(7, ik, itf) = phinewt/(1.5d0*h0mpc**2*(this%Omega_b*O0_BARYON(this)%density_ratio(a)+this%Omega_c*O0_CDM(this)%density_ratio(a))*a**2)
          endif
          trans(10, ik, itf) =  phiweyl/2.d0  !!check in CAMB, transfer_weyl = 10
          trans(11, ik, itf) = vc/kMpc**2
       enddo
    endif
    deallocate(indices)
  end subroutine coop_cosmology_firstorder_camb_DoSourceK


  subroutine coop_cosmology_firstorder_camb_GetTransfer(this, ik, kMpc, num_trans, tauMpc_trans, trans)
    class(coop_cosmology_firstorder)::this
    type(coop_cosmology_firstorder_source)::s
    COOP_INT ik
    COOP_INT::num_trans
    COOP_REAL,dimension(:)::tauMpc_trans
    COOP_SINGLE,dimension(:,:,:)::trans
    
    COOP_REAL kMpc, h0mpc, psi, phinewt, phiweyl, h0tau, vc
    COOP_INT:: i, itf
    COOP_INT,dimension(:),allocatable::indices
    COOP_REAL::  a, kbyH0
    allocate(indices(num_trans))
    h0mpc = this%H0Mpc()
    kbyH0 = kMpc/h0mpc
    call this%allocate_source(m = 0, source = s, k = (/ kbyH0 /), tau = tauMpc_trans*h0mpc, indices=indices)
    call this%compute_source_k(s, 1, transfer_only = .true.)
    trans(1, ik, :) = kMpc/this%h()
    do itf = 1, num_trans
       psi = s%saux(3, 1, indices(itf))
       phiweyl = s%saux(2, 1, indices(itf))
       h0tau = tauMpc_trans(itf)*h0mpc
       vc = s%saux(6, 1, indices(itf))
       a = this%aoftau(h0tau)
       phinewt = phiweyl - psi
       if(this%index_massivenu.ne.0)then
          trans(7, ik, itf) = phinewt/(1.5d0*h0mpc**2*(this%Omega_b*O0_BARYON(this)%density_ratio(a)+this%Omega_c*O0_CDM(this)%density_ratio(a) + this%Omega_massivenu*O0_MASSIVENU(this)%density_ratio(a) )*a**2)
       else
          trans(7, ik, itf) = phinewt/(1.5d0*h0mpc**2*(this%Omega_b*O0_BARYON(this)%density_ratio(a)+this%Omega_c*O0_CDM(this)%density_ratio(a))*a**2)
       endif  !!check in CAMB, transfer_tot
       trans(10, ik, itf) =  phiweyl/2.d0  !!check in CAMB, transfer_weyl = 10
       trans(11, ik, itf) = vc/kMpc**2
    enddo
    deallocate(indices)
  end subroutine coop_cosmology_firstorder_camb_GetTransfer


