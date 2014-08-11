  COOP_INT, optional::genre
  COOP_UNKNOWN_STRING, optional::name
  COOP_INT, optional:: id
  COOP_REAL, optional::Omega
  COOP_REAL, optional::w
  COOP_REAL, optional::cs2
  COOP_REAL, optional::Omega_massless
  type(coop_function),optional::fwp1, fcs2, fwp1eff
  COOP_INT i
  COOP_REAL_ARRAY::lnrat, warr, cs2arr
  COOP_REAL amin, amax, lnamin, lnamax, w1, w2, dlna, lnma, lnm, lnrho, lnp, dlnrho, dlnp
  if(present(name))then
     this%name = name
  else
     this%name = "COOP_SPECIES"
  endif
  if(present(id))then
     this%id = id
  else
     this%id = 0
  endif

  if(present(Omega))then
     this%Omega = Omega
  else
     this%Omega = 1
  endif

  if(present(Omega_massless))then
     this%Omega_massless = max(min(Omega_massless, this%Omega), COOP_REAL_OF(0.))
  else
     this%Omega_massless = this%Omega
  endif

  if(present(genre))then
     this%genre =genre
     if(this%genre .eq. COOP_SPECIES_MASSIVE_FERMION .or. this%genre .eq. COOP_SPECIES_MASSIVE_BOSON)then
        if(this%Omega_massless .ge. this%Omega * 0.9999)then
           this%genre = COOP_SPECIES_MASSLESS
        endif
     endif
  else
     if(present(Omega_massless))then
        if(Omega_massless .le. 0.)then
           this%genre = COOP_SPECIES_CDM
        elseif(Omega_massless/this%Omega .le. 0.9999)then
           this%genre = COOP_SPECIES_MASSIVE_FERMION
        else
           this%genre = COOP_SPECIES_MASSLESS
        endif
     else
        this%genre = COOP_SPECIES_FLUID
     endif
  endif

  select case(this%genre)
  case(COOP_SPECIES_MASSLESS)
     this%wp1 = 1.d0+1.d0/3.
     this%cs2 = 1.d0/3.d0
  case(COOP_SPECIES_CDM)
     this%wp1 = 1.d0
     this%cs2 = 0.d0
  case(COOP_SPECIES_COSMOLOGICAL_CONSTANT)
     this%wp1 = 0.d0
     this%cs2 = 1.d0
  case default
     if(present(w))then
        this%wp1 = w+1.d0
     else
        this%wp1 = 1.d0
     endif
     if(present(cs2))then
        this%cs2 = cs2
     else
        this%cs2 = 0.d0
     endif
  end select
  select case(this%genre)
  case(COOP_SPECIES_MASSLESS, COOP_SPECIES_COSMOLOGICAL_CONSTANT, COOP_SPECIES_CDM)
     this%w_dynamic = .false.
     this%cs2_dynamic = .false.
     return
  case(COOP_SPECIES_MASSIVE_BOSON)
     this%w_dynamic = .true.
     this%cs2_dynamic = .true.
     amin = coop_min_scale_factor
     amax = coop_scale_factor_today
     lnamin = log(amin)
     lnamax = log(amax)
     dlna = (lnamax - lnamin)/(coop_default_array_size  - 1)
     w1 = log(this%Omega/this%Omega_massless)
     call coop_boson_get_lnam(w1, lnm)
     this%mbyT = exp(lnm)
     !$omp parallel do private(lnma, lnrho, lnp, dlnrho, dlnp)
     do i=1, coop_default_array_size
        lnma = lnm + lnamin+(i-1)*dlna
        call coop_boson_get_lnrho(lnma, lnrho)
        call coop_boson_get_lnp(lnma, lnp)
        call coop_boson_get_dlnrho(lnma, dlnrho)
        call coop_boson_get_dlnp(lnma, dlnp)
        warr(i) = exp(lnp-lnrho)
        cs2arr(i) = exp(dlnp-dlnrho)*warr(i)
        lnrat(i) = lnrho- 4.d0*(lnamin+(i-1)*dlna) - w1
     enddo
     !$omp end parallel do
     call this%fwp1%init(coop_default_array_size, amin, amax, 1.d0+warr, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false.)
     call this%fcs2%init(coop_default_array_size, amin, amax, cs2arr, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false.)
     call this%flnrho%init(coop_default_array_size, amin, amax, lnrat, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false.)
     this%fwp1eff = this%fwp1
  case(COOP_SPECIES_MASSIVE_FERMION)
     this%w_dynamic = .true.
     this%cs2_dynamic = .true.
     amin = coop_min_scale_factor
     amax = coop_scale_factor_today
     lnamin = log(amin)
     lnamax = log(amax)
     dlna = (lnamax - lnamin)/(coop_default_array_size  - 1)
     w1 = log(this%Omega/this%Omega_massless)
     call coop_fermion_get_lnam( w1 , lnm)
     this%mbyT = exp(lnm)
     !$omp parallel do private(lnma, lnrho, lnp, dlnrho, dlnp)
     do i=1, coop_default_array_size
        lnma = lnm + lnamin+(i-1)*dlna
        call coop_fermion_get_lnrho(lnma, lnrho)
        call coop_fermion_get_lnp(lnma, lnp)
        call coop_fermion_get_dlnrho(lnma, dlnrho)
        call coop_fermion_get_dlnp(lnma, dlnp)
        warr(i) = exp(lnp-lnrho)
        cs2arr(i) = exp(dlnp-dlnrho)*warr(i)
        lnrat(i) = lnrho- 4.d0*(lnamin+(i-1)*dlna) - w1
     enddo
     !$omp end parallel do
     call this%fwp1%init(coop_default_array_size, amin, amax, 1.d0+warr, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false.)
     call this%fcs2%init(coop_default_array_size, amin, amax, cs2arr, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false.)
     call this%flnrho%init(coop_default_array_size, amin, amax, lnrat, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false.)
     this%fwp1eff = this%fwp1
  case default
     if(present(fwp1) .or. present(fwp1eff))then
        this%w_dynamic = .true.
        if(present(fwp1))then
           this%fwp1 = fwp1
        else
           if(present(w))then              
              warr =  w
              call this%fwp1%init(coop_default_array_size, coop_min_scale_factor, coop_scale_factor_today, 1.d0+warr, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false.)
           else
              stop "For dynamic w you have to pass either w or fwp1 to species_init"
           endif
        endif
        if(present(fwp1eff))then
           this%fwp1eff = fwp1eff
        else
           this%fwp1eff = fwp1
        endif
     else
        this%w_dynamic = .false.
     endif
     if(present(fcs2))then
        this%cs2_dynamic = .true.
        this%fcs2 = fcs2       
     else
        this%cs2_dynamic = .false.
     endif
     if(this%w_dynamic)then
        amin = coop_min_scale_factor
        amax = coop_scale_factor_today
        lnamin = log(amin)
        lnamax = log(amax)
        dlna = (lnamax - lnamin)/(coop_default_array_size  - 1)
        w1 = this%wp1effofa(amax)
        lnrat(coop_default_array_size) = 0.d0
        do i= coop_default_array_size - 1, 1, -1
           w2 = this%wp1effofa(exp((lnamin + dlna*(i-1))))
           lnrat(i) = lnrat(i+1) - (w2 + w1 + 4.*this%wp1effofa(exp((lnamin + dlna*(i-0.5)))))
           w1 = w2
        enddo
        lnrat = lnrat*(-dlna/2.)
        call this%flnrho%init(coop_default_array_size, amin, amax, lnrat, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false.)
     endif
  end select
  if(this%w_dynamic)then
     this%wp1 = this%wp1ofa(coop_scale_factor_today)
     call this%flnrho%set_boundary(slopeleft = -3.*this%wp1effofa(amin), sloperight = -3.*this%wp1effofa(amax))
  endif
  if(this%cs2_dynamic) this%cs2 = this%cs2ofa(coop_scale_factor_today)
