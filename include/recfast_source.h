!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!C Integrator for Cosmic Recombination of Hydrogen and Helium,
!!C developed by Douglas Scott (dscott@astro.ubc.ca)
!!C based on calculations in the papers Seager, Sasselov & Scott
!!C (ApJ, 523, L1, 1999; ApJS, 128, 407, 2000)
!!C and "fudge" updates in Wong, Moss & Scott (2008).
!!C
!!C Permission to use, copy, modify and distribute without fee or royalty at
!!C any tier, this software and its documentation, for any purpose and without
!!C fee or royalty is hereby granted, provided that you agree to comply with
!!C the following copyright notice and statements, including the disclaimer,
!!C and that the same appear on ALL copies of the software and documentation,
!!C including modifications that you make for internal use or for distribution:
!!C
!!C Copyright 1999-2010 by University of British Columbia.  All rights reserved.
!!C
!!C THIS SOFTWARE IS PROVIDED "AS IS", AND U.B.C. MAKES NO 
!!C REPRESENTATcoop_recfast_ionS OR WARRANTIES, EXPRESS OR IMPLIED.  
!!C BY WAY OF EXAMPLE, BUT NOT LIMITATcoop_recfast_ion,
!!C U.B.C. MAKES NO REPRESENTATcoop_recfast_ionS OR WARRANTIES OF 
!!C MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT 
!!C THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATcoop_recfast_ion WILL NOT INFRINGE 
!!C ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.   
!!C
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!CN	Name:	RECFAST
!!CV	Version: 1.5.2
!!C 
!!CP	Purpose:  Calculate ionised fraction as a function of redshift.
!!CP		  Solves for H and He simultaneously, and includes
!!CP		  H "fudge factor" for low z effect, as well as
!!CP	          HeI fudge factor.
!!C
!!CD	Description: Solves for ionisation history since recombination
!!CD	using the equations in Seager, Sasselov & Scott (ApJ, 1999).
!!CD	The Cosmological model can be flat or open.
!!CD	The matter temperature is also followed, with an update from
!!CD	Moss & Scott (2009).
!!CD	The values for \alpha_B for H are from Hummer (1994).
!!CD	The singlet HeI coefficient is a fit from the full code.
!!CD	Additional He "fudge factors" are as described in Wong, Moss
!!CD	and Scott (2008).
!!CD	Extra fitting function included (in optical depth) to account
!!CD	for extra H physics described in Rubino-Martin et al. (2010).
!!CD	Care is taken to use the most accurate constants.
!!CD	Note that some parameters are fixed (e.g. N_nu=3, nu's are
!!CD	massless, w=-1, etc.) - some users may want to explictly
!!CD	imput their own H(z) to account for extra physics.
!!CD	This is provided as a PROGRAM, which can be easily converted
!!CD	to a SUBROUTINE for use in CMB Boltzmann codes.
!!C		
!!CA	Arguments:
!!CA	Name, Description
!!CA	Double precision throughout
!!CA
!!CA	z is redshift - W is sqrt(1+z), like conformal time
!!CA	x is total ionised fraction, relative to H
!!CA	x_H is ionized fraction of H - y(1) in R-K routine
!!CA	x_He is ionized fraction of He - y(2) in R-K routine
!!CA	(note that x_He=n_He+/n_He here and not n_He+/n_H)
!!CA	Tmat is matter temperature - y(3) in R-K routine
!!CA	f's are the derivatives of the Y's
!!CA	alphaB is case B recombination rate
!!CA	alpHe is the singlet only HeII recombination rate
!!CA	a_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
!!CA	b_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
!!CA	c_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
!!CA	d_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
!!CA	a_VF is Verner and Ferland type fitting parameter for Helium
!!CA	b_VF is Verner and Ferland type fitting parameter for Helium
!!CA	T_0 is Verner and Ferland type fitting parameter for Helium
!!CA	T_1 is Verner and Ferland type fitting parameter for Helium
!!CA	Tnow is the observed CMB temperature today
!!CA	OmegaT is the total Omega_0
!!CA      OmegaL is the Omega_0 contribution from a Cosmological constant
!!CA      OmegaK is the Omega_0 contribution in curvature (1-O_T-O_L)
!!CA      OmegaB is Omega in baryons today
!!CA	OmegaC is the Omega_0 in (cold) dark matter: OmegaT=OmegaC+OmegaB
!!CA	Yp is the primordial helium abundace
!!CA	fHe is He/H number ratio = Yp/4(1-Yp)
!!CA	Trad and Tmat are radiation and matter temperatures
!!CA	epsilon is the approximate difference (=Trad-Tmat) at high z
!!CA	OmegaB is Omega in baryons today
!!CA	H is Hubble constant in units of 100 km/s/Mpc
!!CA	HOinp is input value of Hubble constant in units of 100 km/s/Mpc
!!CA	HO is Hubble constant in SI units
!!CA	recfast_bigH is 100 km/s/Mpc in SI units
!!CA	Hz is the value of H at the specific z (in coop_recfast_ion)
!!CA	G is grvitational constant
!!CA	n is number density of hydrogen
!!CA	Nnow is number density today
!!CA	x0 is initial ionized fraction
!!CA	x_H0 is initial ionized fraction of Hydrogen
!!CA	x_He0 is initial ionized fraction of Helium
!!CA	rhs is dummy for calculating x0
!!CA	recfast_z_initial and recfast_z_final are starting and ending redshifts
!!CA	fnu is the contribution of neutrinos to the radn. energy density
!!CA	zeq is the redshift of matter-radiation equality
!!CA	zstart and zend are for each pass to the integrator
!!CA	w0 and w1 are conformal-time-like initial and final zi and zf's
!!CA	Lw0 and Lw1 are logs of w0 and w1
!!CA	hw is the interval in W
!!CA	C,coop_SI_kb,coop_SI_h: speed of light, Boltzmann's and Planck's constants
!!CA	coop_SI_m_e,coop_SI_m_H: electron mass and H atomic mass in SI
!!CA	recfast_not4: ratio of 4He atomic mass to 1H atomic mass
!!CA	sigma: Thomson cross-section
!!CA	a: radiation constant for u=aT^4
!!CA	coop_pi: coop_pi
!!CA	A2s1s: 2s-1s two photon rate for Hydrogen
!!CA	A2s1s_He: 2s-1s two photon rate for Helium
!!CA	DeltaB: energy of first excited state from continuum = 3.4eV
!!CA	DeltaB_He: energy of first excited state from cont. for He = 3.4eV
!!CA	recfast_L_H_ion: level for H ionization in m^-1
!!CA	recfast_L_H_alpha: level for H Ly alpha in m^-1
!!CA	recfast_L_He1_ion: level for HeI ionization
!!CA	recfast_L_He2_ion: level for HeII ionization
!!CA	recfast_L_He_2s: level for HeI 2s
!!CA	recfast_L_He_2p: level for He 2p (21P1-11S0) in m^-1
!!CA	recfast_Lalpha: Ly alpha wavelength in SI
!!CA	recfast_Lalpha_He: Helium I 2p-1s wavelength in SI
!!CA	mu_H,mu_T: mass per H atom and mass per particle
!!CA	H_frac: follow Tmat when t_Compton / t_Hubble > H_frac
!!CA	dHdz is the derivative of H at the specific z (in coop_recfast_ion)
!!CA	recfast_CDB=DeltaB/coop_SI_kb			Constants derived from B1,B2,R
!!CA	recfast_CDB_He=DeltaB_He/coop_SI_kb		n=2-infinity for He in Kelvin
!!CA	recfast_CB1=recfast_CDB*4.			recfast_Lalpha and sigma_Th, calculated
!!CA	recfast_CB1_He1: recfast_CB1 for HeI ionization potential
!!CA	recfast_CB1_He2: recfast_CB1 for HeII ionization potential
!!CA	CR=2*coop_pi*(coop_SI_m_e/coop_SI_h)*(coop_SI_kb/coop_SI_h)	once and passed in a common block
!!CA	recfast_CK=recfast_Lalpha**3/(8.*coop_pi)
!!CA	recfast_CK_He=recfast_Lalpha_He**3/(8.*coop_pi)
!!CA	recfast_CL=C*coop_SI_h/(coop_SI_kb*recfast_Lalpha)
!!CA	recfast_Crecfast_L_He=C*coop_SI_h/(coop_SI_kb*recfast_Lalpha_He)
!!CA	recfast_CT=(8./3.)*(sigma/(coop_SI_m_e*C))*a
!!CA	Bfact=exp((E_2p-E_2s)/kT)	Extra Boltzmann factor
!!CA	fu is a "fudge factor" for H, to approximate low z behaviour
!!CA	b_He is a "fudge factor" for HeI, to approximate higher z behaviour
!!CA	recfast_Heswitch is an integer for modifying HeI recombination
!!CA	Parameters and quantities to describe the extra triplet states
!!CA	 and also the continuum opacity of H, with a fitting function
!!CA	 suggested by KIV, astro-ph/0703438
!!CA	a_trip: used to fit HeI triplet recombination rate
!!CA	b_trip: used to fit HeI triplet recombination rate
!!CA	recfast_L_He_2Pt: level for 23P012-11S0 in m^-1
!!CA	recfast_L_He_2St: level for 23S1-11S0 in m^-1
!!CA	recfast_L_He2St_ion: level for 23S1-continuum in m^-1
!!CA	recfast_A2P_s: Einstein A coefficient for He 21P1-11S0
!!CA	recfast_A2P_t: Einstein A coefficient for He 23P1-11S0    
!!CA	sigma_He_2Ps: H ionization x-section at HeI 21P1-11S0 freq. in m^2
!!CA	sigma_He_2Pt: H ionization x-section at HeI 23P1-11S0 freq. in m^2
!!CA	recfast_CL_PSt = coop_SI_h*coop_SI_c*(recfast_L_He_2Pt - recfast_L_He_2st)/coop_SI_kb
!!CA	CfHe_t: triplet statistical correction
!!CA	recfast_Hswitch is an integer for modifying the H recombination
!!CA	recfast_AGauss1 is the amplitude of the 1st Gaussian for the H fudging
!!CA	recfast_AGauss2 is the amplitude of the 2nd Gaussian for the H fudging
!!CA	recfast_zGauss1 is the ln(1+z) central value of the 1st Gaussian
!!CA	recfast_zGauss2 is the ln(1+z) central value of the 2nd Gaussian
!!CA	recfast_wGauss1 is the width of the 1st Gaussian
!!CA	recfast_wGauss2 is the width of the 2nd Gaussian
!!CA	recfast_tol: recfast_tolerance for the integrator
!!CA	cw(24),w(3,9): work space for DVERK
!!CA	Ndim: number of d.e.'s to solve (integer)
!!CA	Nz: number of output redshitf (integer)
!!CA	I: loop index (integer)
!!CA	ind,nw: work-space for DVERK (integer)
!!C
!!CG	Global data (common blocks) referenced:
!!CG	/zLIST/recfast_z_initial,recfast_z_final,Nz
!!CG	/Cfund/C,coop_SI_kb,coop_SI_h,coop_SI_m_e,coop_SI_m_H,recfast_not4,sigma,a,coop_pi
!!CG	/data/A2s1s,H_frac,recfast_CB1,recfast_CDB,CR,recfast_CK,recfast_CL,recfast_CT,
!!CG		fHe,recfast_CB1_He1,recfast_CB1_He2,recfast_CDB_He,A2s1s_He,Bfact,recfast_CK_He,recfast_Crecfast_L_He
!!CG      /Cosmo/Tnow,HO,Nnow,z_eq,OmegaT,OmegaL,OmegaK
!!CG	/Hemod/b_He,recfast_A2P_s,recfast_A2P_t,sigma_He_2Ps,sigma_He_2Pt,
!!CG		recfast_L_He_2p,recfast_L_He_2Pt,recfast_L_He_2St,recfast_L_He2St_ion
!!CG	/Hmod/recfast_AGauss1,recfast_AGauss2,recfast_zGauss1,recfast_zGauss2,recfast_wGauss1,recfast_wGauss2
!!CG	/Switch/recfast_Heswitch,recfast_Hswitch
!!C
!!CF	File & device access:
!!CF	Unit	/I,IO,O	/Name (if known)
!!C
!!CM	Modules called:
!!CM	DVERK (numerical integrator)
!!CM	GET_INIT (initial values for ionization fractions)
!!CM	coop_recfast_ion (ionization and Temp derivatives)
!!C
!!CC	Comments:
!!CC	none
!!C
!!CH	History:
!!CH	CREATED		(simplest version) 19th March 1989
!!CH	RECREATED	11th January 1995
!!CH			includes variable Cosmology
!!CH			uses DVERK integrator
!!CH			initial conditions are Saha
!!CH	TESTED		a bunch, well, OK, not really
!!CH	MODIFIED	January 1995 (include Hummer's 1994 alpha table)
!!CH			January 1995 (include new value for 2s-1s rate)
!!CH			January 1995 (expand comments)
!!CH			March 1995 (add Saha for Helium)
!!CH			August 1997 (add HeII alpha table)
!!CH			July 1998 (include OmegaT correction and H fudge factor)
!!CH			Nov 1998 (change Trad to Tmat in Rup)
!!CH			Jan 1999 (tidied up for public consumption)
!!CH			Sept 1999 (switch to formula for alpha's, fix glitch)
!!CH			Feb 2000 (fixed overflow problem in He_Boltz)
!!CH			Oct 2001 (fixed OmegaT in z_eq)
!!CH			June 2003 (fixed error in Rdown_He formula)
!!CH			June 2003 (fixed He recombination coefficients)
!!CH			June 2003 (comments to point out fixed N_nu etc.)
!!CH			Oct 2006 (included new value for G)
!!CH			Oct 2006 (improved coop_SI_m_He/coop_SI_m_H to be "recfast_not4")
!!CH			Oct 2006 (fixed error, x for x_H in part of f(1))
!!CH			Jan 2008 (improved HeI recombination effects,
!!CH			              including HeI rec. fudge factor)
!!CH			Feb 2008 (avoid calculating gamma_2Ps and
!!CH			               gamma_2Pt when x_H close to 1.0)
!!CH			Aug 2008 (correction for x_H when Heflag=2
!!CH				     and Helfag>=5 to be smoother)
!!CH			Sept 2008 (added extra term to make transition
!!CH				     smoother for Tmat evolution)
!!CH			Jan 2010 (added fitting function to modify K
!!CH				to match x_e(z) for new H physics)
!!CH			July 2012 (modified fudge factors for better
!!CH				match to codes with more detailed physics)
!!CH			Sept 2012 (fixed "fu" at low z to match modifications)
!!C-
!!C	===============================================================

  Function coop_reionization_xe(z, reionFrac, zre, deltaz)
    COOP_REAL z, coop_reionization_xe, zre, reionFrac, deltaz
    coop_reionization_xe = ReionFrac * (1.+ tanh(((1.+ zre)**1.5 - (1.+z)**1.5)/(1.5*deltaz*sqrt(1.+zre))))/2.d0
  End Function coop_reionization_xe

  
  subroutine coop_recfast_get_xe(bg, xeofa, Tbofa, reionFrac, zre, deltaz)  
    type(coop_cosmology_firstorder)::bg
    type(coop_arguments)::args
    type(coop_function)::xeofa, Tbofa
    integer ndim , nw, ind, i
    logical::reset_space
    integer, parameter::nz = 8192
    COOP_REAL,dimension(nz)::alist,xelist, Tblist, cs2blist
    COOP_REAL reionFrac, zre, deltaz
    !!arguments
    COOP_REAL zend, zstart
    COOP_REAL mu_H, fHe, nnow, HO, fu
    COOP_REAL y(3), rhs, cw(24), w(3, 9), yp(3)

    reset_space = .true.
    ndim = 3
    mu_H = 1.d0/(1.d0-bg%YHe())			!Mass per H atom
    fHe = bg%YHe()/(recfast_not4*(1.d0-bg%YHe()))		!n_He_tot / n_H_tot
    HO = bg%h()*recfast_bigH

    Nnow = 3.d0*(HO**2*coop_Mpsq0) * bg%Omega_b/(8.d0*coop_pi*coop_SI_G*mu_H*coop_SI_m_H)
    !!C	Fudge factor to approximate the low z out of equilibrium effect
    if (recfast_Hswitch .eq. 0) then
       fu=1.14d0
    else
       fu=1.125d0
    end if

    call args%init( r = (/ fHe, Nnow, HO, fu , reionFrac, zre, deltaz /), l = (/ .false. /) )


    call coop_set_uniform(nz, alist, -log(recfast_z_initial), -log(1.d0+recfast_z_final))
    alist = exp(alist)

    !!c	Set up work-space stuff for DVERK
    ind  = 1
    nw   = 3
    cw = 0
    zend = recfast_z_initial
    
    do i = 1, Nz
       zstart = zend
       zend   = 1.d0/alist(i) - 1.d0
       if(args%l(1))then
          call coop_dverk_firstorder_with_args(ndim, coop_recfast_ion, bg, args, zstart,y, zend, recfast_tol, ind, cw, nw, w)
          xelist(i)  = coop_reionization_xe(zend, reionFrac, zre, deltaz)
          Tblist(i) = y(3)
       elseif (zend.gt.8000.d0) then
          xelist(i) = 1.d0+2.d0*fHe
          Tblist(i) = bg%Tcmb()*(1.d0+zend)
       else if(zend.gt.5000.d0)then
          rhs = dexp( 1.5d0 * dLog(recfast_CR*bg%Tcmb()/(1.d0+zend)) & 
               - recfast_CB1_He2/(bg%Tcmb()*(1.d0+zend)) ) / Nnow
          xelist(i) = 0.5d0 * ( dsqrt( (rhs-1.d0-fHe)**2 & 
               + 4.d0*(1.d0+2.d0*fHe)*rhs) - (rhs-1.d0-fHe) )
          Tblist(i) = bg%Tcmb()*(1.d0+zend)
       else if(zend.gt.3500.d0)then
          xelist(i) = 1.d0 + fHe
          y(1) = 1.d0
          y(2) = 1.d0
          y(3) = bg%Tcmb()*(1.d0+zend)
          Tblist(i) = y(3)
       else if(y(2).gt. recfast_He_truncate)then

          rhs = dexp( 1.5d0 * dLog(recfast_CR*bg%Tcmb()/(1.d0+zend)) & 
               - recfast_CB1_He1/(bg%Tcmb()*(1.d0+zend)) ) / Nnow
          rhs = rhs*4.d0		!ratio of g's is 4 for He+ <-> He0
          xelist(i) = 0.5d0 * ( dsqrt( (rhs-1.d0)**2 + 4.d0*(1.d0+fHe)*rhs ) & 
               - (rhs-1.d0))
          y(1) = 1.d0
          y(2) = (xelist(i) - 1.d0)/fHe
          y(3) = bg%Tcmb()*(1.d0+zend)
          Tblist(i) = y(3)
       else if (y(1).gt. recfast_H_truncate) then
          call coop_dverk_firstorder_with_args(ndim, coop_recfast_ion, bg, args, zstart,  y, zend, recfast_tol, ind, cw, nw, w)
          rhs = dexp( 1.5d0 * dLog(recfast_CR*bg%Tcmb()/(1.d0+zend)) & 
               - recfast_CB1/(bg%Tcmb()*(1.d0+zend)) ) / Nnow
          y(1) = 0.5d0 * (dsqrt( rhs**2+4.d0*rhs ) - rhs )
          xelist(i) = y(1) + fHe*y(2)
          Tblist(i) = y(3)
       else
          if(reset_space)then
             ind = 1
             w = 0
             cw = 0
             reset_space = .false.
          endif
          call coop_dverk_firstorder_with_args(ndim, coop_recfast_ion, bg, args, zstart,y, zend, recfast_tol, ind, cw, nw, w)
          xelist(i)  = y(1) + fHe*y(2)
          Tblist(i) = y(3)
          if(zend .gt. zre - deltaz*10.d0)then
             if(xelist(i) .lt.  coop_reionization_xe(zend, reionFrac, zre, deltaz))then
                xelist(i) = coop_reionization_xe(zend, reionFrac, zre, deltaz)
                args%l(1) = .true.
                ind = 1
                cw = 0
                w = 0
             endif
          endif
       end if
    end do

    call xeofa%init(n = nz, xmin=alist(1), xmax = alist(nz), f = xelist, xlog = .true., ylog = .true., fleft = xelist(1), fright = xelist(nz), slopeleft= 0.d0, sloperight = 0.d0, check_boundary = .false., method = COOP_INTERPOLATE_QUADRATIC)

    call Tbofa%init(n = nz, xmin=alist(1), xmax = alist(nz), f = Tblist, xlog = .true., ylog = .true., fleft = Tblist(1), fright = Tblist(nz), slopeleft= -1.d0, sloperight = 0.d0, check_boundary = .false., method = COOP_INTERPOLATE_QUADRATIC)
    do i=1, nz
       cs2blist(i) =max(coop_SI_barssc0*(1.d0-0.75d0*bg%Yhe()+(1.d0-bg%YHe())*xelist(i))  &
            *tblist(i)*(1.d0 - Tbofa%derivative_bare(log(alist(i)))/3.d0), 1.d-20)
    enddo

    call bg%species(bg%index_baryon)%fcs2%init(n = nz, xmin=alist(1), xmax = alist(nz), f = cs2blist, xlog = .true., ylog = .true., fleft = cs2blist(1), fright = cs2blist(nz), slopeleft= -1.d0, sloperight = 0.d0, check_boundary = .false., method = COOP_INTERPOLATE_LINEAR)
    call args%free()
  end subroutine coop_recfast_get_xe

  subroutine coop_recfast_ion(Ndim, z, y, f, bg, args)
    type(coop_cosmology_firstorder)::bg
    type(coop_arguments)::args
    integer Ndim,Heflag
    COOP_REAL, parameter:: a_PPB = 4.309d0, b_PPB = -0.6166d0,     c_PPB = 0.6703d0,   d_PPB = 0.5300d0,  a_VF = 10.d0**(-16.744d0),    b_VF = 0.711d0,    T_0 = 10.d0**(0.477121d0)	,   T_1 = 10.d0**(5.114d0),     a_trip = 10.d0**(-16.306d0),    b_trip = 0.761D0
    COOP_REAL x_H, x_He, x, Tmat, n, n_He, Trad, Hz, dHdz, Rdown, Rup, sq_0, sq_1, Rdown_He, Rup_He, He_Boltz, K, Rdown_trip, Rup_trip, K_He,  tauHe_s, pHe_s, Doppler, AHcon, gamma_2Ps,  CfHe_t,recfast_CL_PSt, epsilon, factor, qb, pb, timeTh, timeH,  gamma_2Pt ,  pHe_t , tauHe_t
    COOP_REAL z, y(3), f(3), scal
    

#define ARGS_FHE args%r(1)
#define ARGS_NNOW args%r(2)
#define ARGS_H0  args%r(3)
#define ARGS_FUDGE args%r(4)
#define ARGS_REION_FRAC args%r(5)
#define ARGS_REION_ZRE args%r(6)
#define ARGS_REION_DELTAZ args%r(7)
#define ARGS_REION_SWITCH args%l(1)

    scal = 1.d0/(1.d0+z)
    Hz = ARGS_H0 * bg%Hratio(scal)
    !!c	Also calculate derivative for use later
    dHdz = -bg%HdotbyHsq(scal)*Hz*scal



    Trad = bg%Tcmb() * (1.d0+z)
    Tmat = y(3)
    if(ARGS_REION_SWITCH)then
       f(1) = 0
       f(2) = 0
       x =  coop_reionization_xe(z, ARGS_REION_FRAC, ARGS_REION_ZRE, ARGS_REION_DELTAZ) 
       f(3)= recfast_CT * (Trad**4) * x / (1.d0+x+ARGS_FHE) & 
            * (Tmat-Trad) / (Hz*(1.d0+z)) + 2.d0*Tmat/(1.d0+z)
       return
    endif

    x_H = y(1)
    x_He = y(2)
    x = x_H + ARGS_FHE * x_He
    n = ARGS_NNOW * bg%species(bg%index_baryon)%density_ratio(scal)
    n_He = ARGS_FHE * n

    !!c	Get the radiative rates using PPQ fit (identical to Hummer's table)
    Rdown=1.d-19*a_PPB*(Tmat/1.d4)**b_PPB & 
         /(1.d0+c_PPB*(Tmat/1.d4)**d_PPB)
    Rup = Rdown * (recfast_CR*Tmat)**(1.5d0)*dexp(-recfast_CDB/Tmat)

    !!c	calculate He using a fit to a Verner & Ferland type formula
    sq_0 = dsqrt(Tmat/T_0)
    sq_1 = dsqrt(Tmat/T_1)
    !!c	typo here corrected by Wayne Hu and Savita Gahlaut
    Rdown_He = a_VF/(sq_0*(1.d0+sq_0)**(1.d0-b_VF))
    Rdown_He = Rdown_He/(1.d0+sq_1)**(1.d0+b_VF)
    Rup_He = Rdown_He*(recfast_CR*Tmat)**(1.5d0)*dexp(-recfast_CDB_He/Tmat)
    Rup_He = 4.d0*Rup_He	!statistical weights factor for HeI
    !!c	Avoid overflow (pointed out by Jacques Roland)
    He_Boltz = dexp(min(recfast_Bfact/Tmat, 680.d0))

    !!c	now deal with H and its fudges
    if (recfast_Hswitch.eq.0) then 
       K = recfast_CK/Hz		!Peebles coefficient K=lambda_a^3/8piH
    else
       !!c	fit a double Gaussian correction function
       K = recfast_CK/Hz*(1.0d0 & 
            +recfast_AGauss1*dexp(-((log(1.0d0+z)-recfast_zGauss1)/recfast_wGauss1)**2.d0) & 
            +recfast_AGauss2*dexp(-((log(1.0d0+z)-recfast_zGauss2)/recfast_wGauss2)**2.d0))
    end if

    !!c	add the HeI part, using same T_0 and T_1 values
    Rdown_trip = a_trip / (sq_0*(1.d0+sq_0)**(1.0-b_trip))
    Rdown_trip = Rdown_trip / ((1.d0+sq_1)**(1.d0+b_trip))
    Rup_trip = Rdown_trip * exp(- coop_SI_h * coop_SI_c*recfast_L_He2St_ion/(coop_SI_kb*Tmat))
    Rup_trip = Rup_trip*((recfast_CR*Tmat)**(1.5d0))*(4.d0/3.d0)
    !!c	last factor here is the statistical weight

    !!c       try to avoid "NaN" when x_He gets too small
    if ((x_He.lt.5.d-9) .or. (x_He.gt.0.980)) then 
       Heflag = 0
    else
       Heflag = recfast_Heswitch
    end if
    if (Heflag.eq.0)then		!use Peebles coeff. for He
       K_He = recfast_CK_He/Hz
    else	!for Heflag>0 		!use Sobolev escape probability
       tauHe_s = recfast_A2P_s*recfast_CK_He*3.d0*n_He*(1.d0-x_He)/Hz
       pHe_s = (1.d0 - dexp(-tauHe_s))/tauHe_s
       K_He = 1.d0/(recfast_A2P_s*pHe_s*3.d0*n_He*(1.d0-x_He))
       !!c	smoother criterion here from Antony Lewis & Chad Fendt
       if (((Heflag.eq.2).or.(Heflag.ge.5)).and.(x_H.lt.0.9999999d0))then
          !!c	use fitting formula for continuum opacity of H
          !!c	first get the Doppler width parameter
          Doppler = 2.D0*coop_SI_kb*Tmat/(coop_SI_m_H*recfast_not4*coop_SI_c*coop_SI_c)
          Doppler = coop_SI_c*recfast_L_He_2p*dsqrt(Doppler)
          gamma_2Ps = 3.d0*recfast_A2P_s*ARGS_FHE*(1.d0-x_He)*coop_SI_c*coop_SI_c & 
               /(dsqrt(coop_pi)*sigma_He_2Ps*8.d0*coop_pi*Doppler*(1.d0-x_H)) & 
               /((coop_SI_c*recfast_L_He_2p)**2.d0)
          pb = 0.36d0  !value from KIV (2007)
          qb = recfast_b_He

          !!c	calculate AHcon, the value of A*p_(con,H) for H continuum opacity
          AHcon = recfast_A2P_s/(1.d0+pb*(gamma_2Ps**qb))
          K_He=1.d0/((recfast_A2P_s*pHe_s+AHcon)*3.d0*n_He*(1.d0-x_He))
       end if
       if (Heflag.ge.3) then		!include triplet effects
          tauHe_t = recfast_A2P_t*n_He*(1.d0-x_He)*3.d0
          tauHe_t = tauHe_t /(8.d0*coop_pi*Hz*recfast_L_He_2Pt**(3.d0))
          pHe_t = (1.d0 - dexp(-tauHe_t))/tauHe_t
          recfast_CL_PSt = coop_SI_h*coop_SI_c*(recfast_L_He_2Pt - recfast_L_He_2st)/coop_SI_kb
	  if ((Heflag.eq.3) .or. (Heflag.eq.5).or.(x_H.gt.0.99999d0)) then
      !!c	no H cont. effect
             CfHe_t = recfast_A2P_t*pHe_t*dexp(-recfast_CL_PSt/Tmat)
             CfHe_t = CfHe_t/(Rup_trip+CfHe_t)	!"C" factor for triplets
	  else					!include H cont. effect
             Doppler = 2.d0*coop_SI_kb*Tmat/(coop_SI_m_H*recfast_not4*coop_SI_c*coop_SI_c)
             Doppler = coop_SI_c*recfast_L_He_2Pt*dsqrt(Doppler)
             gamma_2Pt = 3.d0*recfast_A2P_t*ARGS_FHE*(1.d0-x_He)*coop_SI_c*coop_SI_c & 
                  /(dsqrt(coop_pi)*sigma_He_2Pt*8.d0*coop_pi*Doppler*(1.d0-x_H)) & 
                  /((coop_SI_c*recfast_L_He_2Pt)**2.d0)
             !!c	use the fitting parameters from KIV (2007) in this case
             pb = 0.66d0
             qb = 0.9d0
             AHcon = recfast_A2P_t/(1.d0+pb*gamma_2Pt**qb)/3.d0
             CfHe_t = (recfast_A2P_t*pHe_t+AHcon)*dexp(-recfast_CL_PSt/Tmat)
             CfHe_t = CfHe_t/(Rup_trip+CfHe_t)	!"C" factor for triplets
	  end if
       end if
    end if

    !!c	Estimates of Thomson scattering time and Hubble time
    timeTh=(1.d0/(recfast_CT*Trad**4))*(1.d0+x+ARGS_FHE)/x	!Thomson time
    timeH=2.d0/(3.d0*ARGS_H0*(1.d0+z)**1.5)		!Hubble time

    !!c	calculate the derivatives
    !!c	turn on H only for x_H<0.99, and use Saha derivative for 0.98<x_H<0.99
    !!c	(clunky, but seems to work)
    if (x_H.gt. recfast_H_truncate) then			!don't change at all
       f(1) = 0.d0
       !!cc	else if ((x_H.gt.0.98d0).and.(Heflag.eq.0)) then	!don't modify
    else if (x_H.gt.0.985d0) then		!use Saha rate for Hydrogen
       f(1) = (x*x_H*n*Rdown - Rup*(1.d0-x_H)*dexp(-recfast_CL/Tmat)) & 
            /(Hz*(1.d0+z))
       !!c	for interest, calculate the correction factor compared to Saha
       !!c	(without the fudge)
       factor=(1.d0 + K*Recfast_A2s1s*n*(1.d0-x_H)) & 
            /(Hz*(1.d0+z)*(1.d0+K*Recfast_A2s1s*n*(1.d0-x) & 
            +K*Rup*n*(1.d0-x)))
    else					!use full rate for H
       f(1) = ((x*x_H*n*Rdown - Rup*(1.d0-x_H)*dexp(-recfast_CL/Tmat)) & 
            *(1.d0 + K*Recfast_A2s1s*n*(1.d0-x_H))) & 
            /(Hz*(1.d0+z)*(1.d0/ARGS_FUDGE+K*Recfast_A2s1s*n*(1.d0-x_H)/ARGS_FUDGE & 
            +K*Rup*n*(1.d0-x_H)))
    end if
    !!c	turn off the He once it is small
    if (x_He.lt.1.d-15) then
       f(2)=0.d0
    else
       f(2) = ((x*x_He*n*Rdown_He  & 
            - Rup_He*(1.d0-x_He)*dexp(-recfast_Crecfast_L_He/Tmat)) & 
            *(1.d0+ K_He*Recfast_A2s1s_He*n_He*(1.d0-x_He)*He_Boltz)) & 
            /(Hz*(1.d0+z) & 
            * (1.d0 + K_He*(Recfast_A2s1s_He+Rup_He)*n_He*(1.d0-x_He)*He_Boltz))
       !!c	Modification to HeI recombination including channel via triplets
       if (Heflag.ge.3) then
          f(2) = f(2)+ (x*x_He*n*Rdown_trip & 
               - (1.d0-x_He)*3.d0*Rup_trip*dexp(-coop_SI_h*coop_SI_c*recfast_L_He_2st/(coop_SI_kb*Tmat))) & 
               *CfHe_t/(Hz*(1.d0+z))
       end if
    end if

    !!c	follow the matter temperature once it has a chance of diverging

    if (timeTh .lt. recfast_H_frac*timeH) then
       !!c		f(3)=Tmat/(1.d0+z)	!Tmat follows Trad
       !!c	additional term to smooth transition to Tmat evolution,
       !!c	(suggested by Adam Moss)
       epsilon = Hz*(1.d0+x+ARGS_FHE)/(recfast_CT*Trad**3*x)
       f(3) = bg%Tcmb() & 
            + epsilon*((1.d0+ARGS_FHE)/(1.d0+ARGS_FHE+x))*((f(1)+ARGS_FHE*f(2))/x) & 
            - epsilon* dHdz/Hz + 3.0d0*epsilon/(1.d0+z) 
    else
       f(3)= recfast_CT * (Trad**4) * x / (1.d0+x+ARGS_FHE) & 
            * (Tmat-Trad) / (Hz*(1.d0+z)) + 2.d0*Tmat/(1.d0+z)
    end if

    return

  end subroutine coop_recfast_ion



  subroutine coop_cosmology_firstorder_set_zre_from_optre(this)
    class(coop_cosmology_firstorder)::this
    COOP_REAL zremin, zremax, zremid
    COOP_REAL optremin, optremax, optremid, optre_wanted
    integer iloop
    if(this%ReionFrac .le. 0.d0)then
       this%zre = -1.d0
    endif
    optre_wanted = this%optre

    zremin = 0.d0
    this%zre = zremin
    call this%set_optre_from_zre()
    optremin = this%optre

    zremax = 20.d0
    this%zre = zremax
    call this%set_optre_from_zre()    
    optremax = this%optre

    do while(optre_wanted .gt. optremax)
       zremin = zremax
       optremin = optremax
       zremax = zremax * 1.3d0
       this%zre = zremax
       call this%set_optre_from_zre()    
       optremax = this%optre
       if(zremax .gt. 800.d0)then
          call coop_return_error("recfast_set_zre_from_optre", "optical depth too big", "stop")
       endif
    enddo

    iloop = 0
    do while((optremax - optremin) .gt. 1.d-5 .and. iloop .lt. 25)
       zremid = (zremin + zremax)/2.d0
       this%zre = zremid
       call this%set_optre_from_zre()
       optremid = this%optre
       if(optremid .gt. optre_wanted)then
          zremax = zremid
          optremax = optremid
       else
          zremin = zremid
          optremin = optremid
       endif
       iloop = iloop + 1
    enddo
    this%optre = optre_wanted
  end subroutine coop_cosmology_firstorder_set_zre_from_optre

  subroutine coop_cosmology_firstorder_set_optre_from_zre(this)
    class(coop_cosmology_firstorder)::this
    COOP_REAL :: optre, ReionFrac, zre, delta
    if(this%ReionFrac .le. 0.d0)then
       this%optre = 0.d0
    else
       select type(this)
       type is(coop_cosmology_firstorder)
          this%optre = coop_integrate_firstorder(coop_optre_int, 0.d0, this%zre + this%deltaz * 10.d0, this, COOP_REAL_OF(1.e-7))
       class default
          stop "for compatibility with older versions of gfortran optre_from_zre only works for type coop_cosmology_firstorder"
       end select
    endif
  end subroutine coop_cosmology_firstorder_set_optre_from_zre

  function coop_optre_int(z, cosmology) result(dkappadz)
    COOP_REAL z, dkappadz
    type(coop_cosmology_firstorder) cosmology
    dkappadz = cosmology%dkappadtau_coef * coop_reionization_xe(z, cosmology%reionFrac, cosmology%zre, cosmology%deltaz) / cosmology%Hasq(1.d0/(1.d0+z)) 
  end function coop_optre_int


  subroutine coop_dverk_firstorder_with_args(n, fcn, cosmology, args, x, y, xend, tol, ind, c, nw, w)
    type(coop_cosmology_firstorder) cosmology
    type(coop_arguments) args
#define DVERK_ARGUMENTS ,cosmology,args
#include "dverk_nostop.h"    
#undef DVERK_ARGUMENTS
  end subroutine coop_dverk_firstorder_with_args


  subroutine coop_dverk_firstorder(n, fcn, cosmology, pert, x, y, xend, tol, ind, c, nw, w)
    type(coop_cosmology_firstorder) cosmology
    type(coop_pert_object) pert
#define DVERK_ARGUMENTS ,cosmology,pert
#include "dverk_nostop.h"    
#undef DVERK_ARGUMENTS
  end subroutine coop_dverk_firstorder


  function coop_integrate_firstorder(func, a, b, cosmology, precision) result(integral)
    type(coop_cosmology_firstorder) cosmology
#define QROMB_ARGUMENTS ,cosmology
#include "qromb.h"
#undef QROMB_ARGUMENTS
  end function coop_integrate_firstorder



