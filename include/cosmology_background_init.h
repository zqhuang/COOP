    this%num_species = 0
    this%Omega_k_value = 1.
    this%need_setup_background = .true.
    if(present(name))then
       this%name = name
    else
       this%name = "COOP_COSMOLOGY_BACKGROUND"
    endif
    if(present(id))then
       this%id =  id
    else
       this%id = 1
    endif
    if(present(h))then
       this%h_value = h
    else
       this%h_value = COOP_DEFAULT_HUBBLE
    endif
    if(present(Tcmb))Then
       this%Tcmb_value = Tcmb
    else
       this%Tcmb_value = COOP_DEFAULT_TCMB
    endif
    if(present(YHe))Then
       this%YHe_value = YHe
    else
       this%YHe_value = COOP_DEFAULT_YHE
    endif
    if(present(Nnu))then
       this%Nnu_value = Nnu
    else
       this%Nnu_value = COOP_DEFAULT_NUM_NEUTRINO_SPECIES
    endif
