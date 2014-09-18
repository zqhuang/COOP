  COOP_REAL, parameter::coop_recfast_a_initial = 1.d-4
  COOP_REAL,parameter::recfast_z_initial = 1.d0/coop_recfast_a_initial + 1.d0
  integer,parameter::recfast_Hswitch = 1, recfast_Heswitch=6
  COOP_REAL,parameter::recfast_z_final = 0.d0
  !!C	--- Parameter statements
  COOP_REAL, parameter::recfast_bigH=100.0D3/(1.0D6*3.0856775807D16)	!Ho in s-1
  COOP_REAL, parameter::recfast_tol=1.d-5				!Tolerance for R-K
  !!C	--- Data
  !!c	note: neglecting deuterium, making an O(e-5) effect
  COOP_REAL, parameter:: recfast_not4 = 3.9715D0
  COOP_REAL, parameter:: recfast_a = 7.565914D-16
  COOP_REAL, parameter:: A2s1s  = 8.2245809d0
  COOP_REAL, parameter:: A2s1s_He = 51.3d0	!new value from Dalgarno
  COOP_REAL, parameter:: L_H_ion = 1.096787737D7 	!level for H ion. (in m^-1)
  COOP_REAL, parameter:: L_H_alpha	 = 8.225916453D6 !averaged over 2 levels
  COOP_REAL, parameter:: Lalpha = 1.d0/L_H_alpha
  COOP_REAL, parameter:: L_He1_ion	= 1.98310772D7	!from Drake (1993)
  COOP_REAL, parameter:: L_He2_ion	= 4.389088863D7	!from JPhysChemRefData (1987)
  COOP_REAL, parameter:: L_He_2s = 1.66277434D7	!from Drake (1993)
  COOP_REAL, parameter:: L_He_2p = 1.71134891D7	!from Drake (1993)
  !!C	2 photon rates and atomic levels in SI units

  COOP_REAL, parameter:: A2P_s = 1.798287D9    !Morton, Wu & Drake (2006)
  COOP_REAL, parameter:: A2P_t = 177.58D0      !Lach & Pachuski (2001)
  COOP_REAL, parameter:: L_He_2Pt	= 1.690871466D7 !Drake & Morton (2007)
  COOP_REAL, parameter:: L_He_2St = 1.5985597526D7 !Drake & Morton (2007)
  COOP_REAL, parameter:: L_He2St_ion = 3.8454693845D6 !Drake & Morton (2007)
  COOP_REAL, parameter:: sigma_He_2Ps = 1.436289D-22  !Hummer & Storey (1998)
  COOP_REAL, parameter:: sigma_He_2Pt = 1.484872D-22  !Hummer & Storey (1998)
  !!C	Atomic data for HeI 

  COOP_REAL, parameter:: AGauss1 = -0.14D0	!Amplitude of 1st Gaussian
  COOP_REAL, parameter:: AGauss2 = 0.079D0	!Amplitude of 2nd Gaussian
  COOP_REAL, parameter:: zGauss1 = 7.28D0	!ln(1+z) of 1st Gaussian
  COOP_REAL, parameter:: zGauss2 = 6.73D0	!ln(1+z) of 2nd Gaussian
  COOP_REAL, parameter:: wGauss1 = 0.18D0	!Width of 1st Gaussian
  COOP_REAL, parameter:: wGauss2 = 0.33D0	!Width of 2nd Gaussian
  !!C	Gaussian fits for extra H physics (fit by Adam Moss, modified by
  !!C	Antony Lewis)

  COOP_REAL, parameter:: Lalpha_He = 1.d0/L_He_2p
  COOP_REAL, parameter:: recfast_DeltaB = coop_SI_h*coop_SI_c*(L_H_ion-L_H_alpha)
  COOP_REAL, parameter:: CDB = recfast_DeltaB/coop_SI_kB
  COOP_REAL, parameter:: recfast_DeltaB_He = coop_SI_h*coop_SI_c*(L_He1_ion-L_He_2s)	!2s, not 2p
  COOP_REAL, parameter:: CDB_He = recfast_DeltaB_He/coop_SI_kB
  COOP_REAL, parameter:: CB1 = coop_SI_h*coop_SI_c*L_H_ion/coop_SI_kB
  COOP_REAL, parameter:: CB1_He1 = coop_SI_h*coop_SI_c*L_He1_ion/coop_SI_kB	!ionization for HeI
  COOP_REAL, parameter:: CB1_He2 = coop_SI_h*coop_SI_c*L_He2_ion/coop_SI_kB	!ionization for HeII
  COOP_REAL, parameter:: recfast_CR = 2.d0*coop_pi*(coop_SI_m_e/coop_SI_h)*(coop_SI_kB/coop_SI_h)
  COOP_REAL, parameter:: recfast_CK = Lalpha**3/(8.d0*coop_pi)
  COOP_REAL, parameter:: recfast_CK_He = Lalpha_He**3/(8.d0*coop_pi)
  COOP_REAL, parameter:: recfast_CL = coop_SI_c*coop_SI_h/(coop_SI_kB*Lalpha)
  COOP_REAL, parameter:: recfast_CL_He = coop_SI_c*coop_SI_h/(coop_SI_kB/L_He_2s)	!comes from det.bal. of 2s-1s
  COOP_REAL, parameter:: recfast_CT = (8.d0/3.d0)*(coop_SI_sigma_thomson /(coop_SI_m_e*coop_SI_c))*recfast_a
  COOP_REAL, parameter:: Bfact = coop_SI_h*coop_SI_c*(L_He_2p-L_He_2s)/coop_SI_kB
  !!C	Matter departs from radiation when t(Th) > H_frac * t(H)
  !!C	choose some safely small number
  COOP_REAL, parameter::  H_frac = 1.D-3
  COOP_REAL, parameter::    b_He = 0.86
