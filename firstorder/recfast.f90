
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
!!CA	bigH is 100 km/s/Mpc in SI units
!!CA	Hz is the value of H at the specific z (in coop_recfast_ion)
!!CA	G is grvitational constant
!!CA	n is number density of hydrogen
!!CA	Nnow is number density today
!!CA	x0 is initial ionized fraction
!!CA	x_H0 is initial ionized fraction of Hydrogen
!!CA	x_He0 is initial ionized fraction of Helium
!!CA	rhs is dummy for calculating x0
!!CA	zinitial and zfinal are starting and ending redshifts
!!CA	fnu is the contribution of neutrinos to the radn. energy density
!!CA	zeq is the redshift of matter-radiation equality
!!CA	zstart and zend are for each pass to the integrator
!!CA	w0 and w1 are conformal-time-like initial and final zi and zf's
!!CA	Lw0 and Lw1 are logs of w0 and w1
!!CA	hw is the interval in W
!!CA	C,k_B,h_P: speed of light, Boltzmann's and Planck's constants
!!CA	m_e,m_H: electron mass and H atomic mass in SI
!!CA	not4: ratio of 4He atomic mass to 1H atomic mass
!!CA	sigma: Thomson cross-section
!!CA	a: radiation constant for u=aT^4
!!CA	Pi: Pi
!!CA	Lambda: 2s-1s two photon rate for Hydrogen
!!CA	Lambda_He: 2s-1s two photon rate for Helium
!!CA	DeltaB: energy of first excited state from continuum = 3.4eV
!!CA	DeltaB_He: energy of first excited state from cont. for He = 3.4eV
!!CA	L_H_ion: level for H ionization in m^-1
!!CA	L_H_alpha: level for H Ly alpha in m^-1
!!CA	L_He1_ion: level for HeI ionization
!!CA	L_He2_ion: level for HeII ionization
!!CA	L_He_2s: level for HeI 2s
!!CA	L_He_2p: level for He 2p (21P1-11S0) in m^-1
!!CA	Lalpha: Ly alpha wavelength in SI
!!CA	Lalpha_He: Helium I 2p-1s wavelength in SI
!!CA	mu_H,mu_T: mass per H atom and mass per particle
!!CA	H_frac: follow Tmat when t_Compton / t_Hubble > H_frac
!!CA	dHdz is the derivative of H at the specific z (in coop_recfast_ion)
!!CA	CDB=DeltaB/k_B			Constants derived from B1,B2,R
!!CA	CDB_He=DeltaB_He/k_B		n=2-infinity for He in Kelvin
!!CA	CB1=CDB*4.			Lalpha and sigma_Th, calculated
!!CA	CB1_He1: CB1 for HeI ionization potential
!!CA	CB1_He2: CB1 for HeII ionization potential
!!CA	CR=2*Pi*(m_e/h_P)*(k_B/h_P)	once and passed in a common block
!!CA	CK=Lalpha**3/(8.*Pi)
!!CA	CK_He=Lalpha_He**3/(8.*Pi)
!!CA	CL=C*h_P/(k_B*Lalpha)
!!CA	CL_He=C*h_P/(k_B*Lalpha_He)
!!CA	CT=(8./3.)*(sigma/(m_e*C))*a
!!CA	Bfact=exp((E_2p-E_2s)/kT)	Extra Boltzmann factor
!!CA	fu is a "fudge factor" for H, to approximate low z behaviour
!!CA	b_He is a "fudge factor" for HeI, to approximate higher z behaviour
!!CA	Heswitch is an integer for modifying HeI recombination
!!CA	Parameters and quantities to describe the extra triplet states
!!CA	 and also the continuum opacity of H, with a fitting function
!!CA	 suggested by KIV, astro-ph/0703438
!!CA	a_trip: used to fit HeI triplet recombination rate
!!CA	b_trip: used to fit HeI triplet recombination rate
!!CA	L_He_2Pt: level for 23P012-11S0 in m^-1
!!CA	L_He_2St: level for 23S1-11S0 in m^-1
!!CA	L_He2St_ion: level for 23S1-continuum in m^-1
!!CA	A2P_s: Einstein A coefficient for He 21P1-11S0
!!CA	A2P_t: Einstein A coefficient for He 23P1-11S0    
!!CA	sigma_He_2Ps: H ionization x-section at HeI 21P1-11S0 freq. in m^2
!!CA	sigma_He_2Pt: H ionization x-section at HeI 23P1-11S0 freq. in m^2
!!CA	CL_PSt = h_P*C*(L_He_2Pt - L_He_2st)/k_B
!!CA	CfHe_t: triplet statistical correction
!!CA	Hswitch is an integer for modifying the H recombination
!!CA	AGauss1 is the amplitude of the 1st Gaussian for the H fudging
!!CA	AGauss2 is the amplitude of the 2nd Gaussian for the H fudging
!!CA	zGauss1 is the ln(1+z) central value of the 1st Gaussian
!!CA	zGauss2 is the ln(1+z) central value of the 2nd Gaussian
!!CA	wGauss1 is the width of the 1st Gaussian
!!CA	wGauss2 is the width of the 2nd Gaussian
!!CA	tol: tolerance for the integrator
!!CA	cw(24),w(3,9): work space for DVERK
!!CA	Ndim: number of d.e.'s to solve (integer)
!!CA	Nz: number of output redshitf (integer)
!!CA	I: loop index (integer)
!!CA	ind,nw: work-space for DVERK (integer)
!!C
!!CG	Global data (common blocks) referenced:
!!CG	/zLIST/zinitial,zfinal,Nz
!!CG	/Cfund/C,k_B,h_P,m_e,m_H,not4,sigma,a,Pi
!!CG	/data/Lambda,H_frac,CB1,CDB,CR,CK,CL,CT,
!!CG		fHe,CB1_He1,CB1_He2,CDB_He,Lambda_He,Bfact,CK_He,CL_He
!!CG      /Cosmo/Tnow,HO,Nnow,z_eq,OmegaT,OmegaL,OmegaK
!!CG	/Hemod/b_He,A2P_s,A2P_t,sigma_He_2Ps,sigma_He_2Pt,
!!CG		L_He_2p,L_He_2Pt,L_He_2St,L_He2St_ion
!!CG	/Hmod/AGauss1,AGauss2,zGauss1,zGauss2,wGauss1,wGauss2
!!CG	/Switch/Heswitch,Hswitch
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
!!CH			Oct 2006 (improved m_He/m_H to be "not4")
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
module coop_recfast_mod
  use coop_wrapper_background
  implicit none
#include "constants.h"

  private


  public:: coop_recfast_get_xe

  integer,parameter::  Hswitch = 1, Heswitch=6

  real*8,parameter::zinitial = 10000.d0
  real*8,parameter::zfinal = 0.d0
  !!C	--- Parameter statements
  real*8, parameter::bigH=100.0D3/(1.0D6*3.0856775807D16)	!Ho in s-1
  real*8, parameter::tol=1.D-5				!Tolerance for R-K


  !!C	--- Data
  real*8, parameter:: C = 2.99792458D8, k_B = 1.380658D-23, h_P = 6.6260755D-34
  real*8, parameter:: m_e = 9.1093897D-31, m_H = 1.673575D-27	!av. H atom
  !!c	note: neglecting deuterium, making an O(e-5) effect
  real*8, parameter:: not4 = 3.9715D0
  real*8, parameter:: sigma = 6.6524616D-29, a = 7.565914D-16
  real*8, parameter:: Pi = 3.141592653589d0
  real*8, parameter:: G = 6.6742D-11
  real*8, parameter:: Lambda  = 8.2245809d0
  real*8, parameter:: Lambda_He	= 51.3d0	!new value from Dalgarno
  real*8, parameter:: L_H_ion = 1.096787737D7 	!level for H ion. (in m^-1)
  real*8, parameter:: L_H_alpha	 = 8.225916453D6 !averaged over 2 levels
  real*8, parameter:: Lalpha = 1.d0/L_H_alpha
  real*8, parameter:: L_He1_ion	= 1.98310772D7	!from Drake (1993)
  real*8, parameter:: L_He2_ion	= 4.389088863D7	!from JPhysChemRefData (1987)
  real*8, parameter:: L_He_2s = 1.66277434D7	!from Drake (1993)
  real*8, parameter:: L_He_2p = 1.71134891D7	!from Drake (1993)
  !!C	2 photon rates and atomic levels in SI units

  real*8, parameter:: A2P_s = 1.798287D9    !Morton, Wu & Drake (2006)
  real*8, parameter:: A2P_t = 177.58D0      !Lach & Pachuski (2001)
  real*8, parameter:: L_He_2Pt	= 1.690871466D7 !Drake & Morton (2007)
  real*8, parameter:: L_He_2St = 1.5985597526D7 !Drake & Morton (2007)
  real*8, parameter:: L_He2St_ion = 3.8454693845D6 !Drake & Morton (2007)
  real*8, parameter:: sigma_He_2Ps = 1.436289D-22  !Hummer & Storey (1998)
  real*8, parameter:: sigma_He_2Pt = 1.484872D-22  !Hummer & Storey (1998)
  !!C	Atomic data for HeI 

  real*8, parameter:: AGauss1 = -0.14D0	!Amplitude of 1st Gaussian
  real*8, parameter:: AGauss2 = 0.079D0	!Amplitude of 2nd Gaussian
  real*8, parameter:: zGauss1 = 7.28D0	!ln(1+z) of 1st Gaussian
  real*8, parameter:: zGauss2 = 6.73D0	!ln(1+z) of 2nd Gaussian
  real*8, parameter:: wGauss1 = 0.18D0	!Width of 1st Gaussian
  real*8, parameter:: wGauss2 = 0.33D0	!Width of 2nd Gaussian
  !!C	Gaussian fits for extra H physics (fit by Adam Moss, modified by
  !!C	Antony Lewis)

  real*8, parameter:: Lalpha_He = 1.d0/L_He_2p
  real*8, parameter:: DeltaB = h_P*C*(L_H_ion-L_H_alpha)
  real*8, parameter:: CDB = DeltaB/k_B
  real*8, parameter:: DeltaB_He = h_P*C*(L_He1_ion-L_He_2s)	!2s, not 2p
  real*8, parameter:: CDB_He = DeltaB_He/k_B
  real*8, parameter:: CB1 = h_P*C*L_H_ion/k_B
  real*8, parameter:: CB1_He1 = h_P*C*L_He1_ion/k_B	!ionization for HeI
  real*8, parameter:: CB1_He2 = h_P*C*L_He2_ion/k_B	!ionization for HeII
  real*8, parameter:: CR = 2.d0*Pi*(m_e/h_P)*(k_B/h_P)
  real*8, parameter:: CK = Lalpha**3/(8.d0*Pi)
  real*8, parameter:: CK_He = Lalpha_He**3/(8.d0*Pi)
  real*8, parameter:: CL = C*h_P/(k_B*Lalpha)
  real*8, parameter:: CL_He = C*h_P/(k_B/L_He_2s)	!comes from det.bal. of 2s-1s
  real*8, parameter:: CT = (8.d0/3.d0)*(sigma/(m_e*C))*a
  real*8, parameter:: Bfact = h_P*C*(L_He_2p-L_He_2s)/k_B
  !!C	Matter departs from radiation when t(Th) > H_frac * t(H)
  !!C	choose some safely small number
  real*8, parameter::  H_frac = 1.D-3
  real*8, parameter::    b_He = 0.86



contains
  
  subroutine coop_recfast_get_xe(bg, xeofa)  
    class(coop_cosmology_background)::bg
    type(coop_arguments)::args
    type(coop_function)::xeofa
    integer ndim , nw, ind, i
    integer, parameter::nz=1000
    COOP_REAL,dimension(nz)::alist,xelist
    !!arguments
    integer index_baryon, index_cdm
    real*8 zend, zstart
    real*8 mu_H, fHe, nnow, HO, fu
    real*8 y(3), rhs, cw(24), w(3, 9), yp(3)
    ndim = 3
    index_baryon = bg%index_of("Baryon")
    index_cdm = bg%index_of("CDM")
    if(index_baryon .eq. 0) stop "Recfast Error: No baryon in the input cosmology"
    if(index_cdm .eq. 0) stop "Recfast Error: No CDM in the input cosmology"
    
    mu_H = 1.d0/(1.d0-bg%YHe())			!Mass per H atom
    fHe = bg%YHe()/(not4*(1.d0-bg%YHe()))		!n_He_tot / n_H_tot
    HO = bg%h()*bigH
    Nnow = 3.d0*(HO)**2*bg%species(index_baryon)%Omega/(8.d0*Pi*G*mu_H*m_H)


    !!C	Fudge factor to approximate the low z out of equilibrium effect
    if (Hswitch .eq. 0) then
       fu=1.14d0
    else
       fu=1.125d0
    end if

    call args%init( i = (/ index_baryon, index_cdm /), r = (/ fHe, Nnow, HO, fu /))


    call coop_set_uniform(nz, alist, -log(zinitial), -log(1.d0+zfinal))
    alist = exp(alist)

    !!c	Set up work-space stuff for DVERK
    ind  = 1
    nw   = 3
    cw = 0
    zend = zinitial
    
    do i = 1, Nz
!!$       zstart = zinitial + (i-1)*(zfinal - zinitial)/nz
!!$       zend =  zinitial + i*(zfinal - zinitial)/nz
       zstart = zend
       zend   = 1.d0/alist(i) - 1.d0
       if (zend.gt.8000.d0) then
          xelist(i) = 1.d0+2.d0*fHe
       else if(zend.gt.5000.d0)then
          rhs = dexp( 1.5d0 * dLog(CR*bg%Tcmb()/(1.d0+zend)) & 
               - CB1_He2/(bg%Tcmb()*(1.d0+zend)) ) / Nnow
          xelist(i) = 0.5d0 * ( dsqrt( (rhs-1.d0-fHe)**2 & 
               + 4.d0*(1.d0+2.d0*fHe)*rhs) - (rhs-1.d0-fHe) )
       else if(zend.gt.3500.d0)then
          xelist(i) = 1.d0 + fHe
          y(1) = 1.d0
          y(2) = 1.d0
          y(3) = bg%Tcmb()*(1.d0+zend)

       else if(y(2).gt.0.99)then

          rhs = dexp( 1.5d0 * dLog(CR*bg%Tcmb()/(1.d0+zend)) & 
               - CB1_He1/(bg%Tcmb()*(1.d0+zend)) ) / Nnow
          rhs = rhs*4.d0		!ratio of g's is 4 for He+ <-> He0
          xelist(i) = 0.5d0 * ( dsqrt( (rhs-1.d0)**2 + 4.d0*(1.d0+fHe)*rhs ) & 
               - (rhs-1.d0))
          y(1) = 1.d0
          y(2) = (xelist(i) - 1.d0)/fHe
          y(3) = bg%Tcmb()*(1.d0+zend)

       else if (y(1).gt.0.99d0) then
          call coop_dverk_with_cosmology(ndim, coop_recfast_ion, bg, args, zstart,  y, zend, tol, ind, cw, nw, w)
          rhs = dexp( 1.5d0 * dLog(CR*bg%Tcmb()/(1.d0+zend)) & 
               - CB1/(bg%Tcmb()*(1.d0+zend)) ) / Nnow
          y(1) = 0.5d0 * (dsqrt( rhs**2+4.d0*rhs ) - rhs )

          xelist(i) = y(1) + fHe*y(2)
       else
          call coop_dverk_with_cosmology(ndim, coop_recfast_ion, bg, args, zstart,y, zend, tol, ind, cw, nw, w)

          xelist(i)  = y(1) + fHe*y(2)
       end if
    end do
    call xeofa%init(n = nz, xmin=alist(1), xmax = alist(nz), f = xelist, xlog = .true., ylog = .false., fleft = xelist(1), fright = xelist(nz), slopeleft= 0.d0, sloperight = 0.d0, check_boundary = .false.)

  end subroutine coop_recfast_get_xe

  subroutine coop_recfast_ion(Ndim, z, y, f, bg, args)
    implicit none
    class(coop_cosmology_background)::bg
    type(coop_arguments)::args
    integer Ndim,Heflag
    real*8, parameter:: a_PPB = 4.309d0, b_PPB = -0.6166d0,     c_PPB = 0.6703d0,   d_PPB = 0.5300d0,  a_VF = 10.d0**(-16.744d0),    b_VF = 0.711d0,    T_0 = 10.d0**(0.477121d0)	,   T_1 = 10.d0**(5.114d0),     a_trip = 10.d0**(-16.306d0),    b_trip = 0.761D0
    real*8 x_H, x_He, x, Tmat, n, n_He, Trad, Hz, dHdz, Rdown, Rup, sq_0, sq_1, Rdown_He, Rup_He, He_Boltz, K, Rdown_trip, Rup_trip, K_He,  tauHe_s, pHe_s, Doppler, AHcon, gamma_2Ps,  CfHe_t,CL_PSt, epsilon, factor, qb, pb, timeTh, timeH,  gamma_2Pt ,  pHe_t , tauHe_t
    real*8 z, y(3), f(3), scal
    

#define ARGS_FHE args%r(1)
#define ARGS_NNOW args%r(2)
#define ARGS_H0  args%r(3)
#define ARGS_FUDGE args%r(4)


    x_H = y(1)
    x_He = y(2)
    x = x_H + ARGS_FHE * x_He
    Tmat = y(3)
    n = ARGS_NNOW * (1.d0+z)**3
    n_He = ARGS_FHE * ARGS_NNOW * (1.d0+z)**3
    Trad = bg%Tcmb() * (1.d0+z)
    scal = 1.d0/(1.d0+z)
    Hz = ARGS_H0 * bg%Hratio(scal)
    !!c	Also calculate derivative for use later
    dHdz = -bg%HdotbyHsq(scal)*Hz/(1.d0+z)
    !!c	Get the radiative rates using PPQ fit (identical to Hummer's table)
    Rdown=1.d-19*a_PPB*(Tmat/1.d4)**b_PPB & 
         /(1.d0+c_PPB*(Tmat/1.d4)**d_PPB)
    Rup = Rdown * (CR*Tmat)**(1.5d0)*dexp(-CDB/Tmat)

    !!c	calculate He using a fit to a Verner & Ferland type formula
    sq_0 = dsqrt(Tmat/T_0)
    sq_1 = dsqrt(Tmat/T_1)
    !!c	typo here corrected by Wayne Hu and Savita Gahlaut
    Rdown_He = a_VF/(sq_0*(1.d0+sq_0)**(1.d0-b_VF))
    Rdown_He = Rdown_He/(1.d0+sq_1)**(1.d0+b_VF)
    Rup_He = Rdown_He*(CR*Tmat)**(1.5d0)*dexp(-CDB_He/Tmat)
    Rup_He = 4.d0*Rup_He	!statistical weights factor for HeI
    !!c	Avoid overflow (pointed out by Jacques Roland)
    He_Boltz = dexp(min(Bfact/Tmat, 680.d0))

    !!c	now deal with H and its fudges
    if (Hswitch.eq.0) then 
       K = CK/Hz		!Peebles coefficient K=lambda_a^3/8piH
    else
       !!c	fit a double Gaussian correction function
       K = CK/Hz*(1.0d0 & 
            +AGauss1*dexp(-((log(1.0d0+z)-zGauss1)/wGauss1)**2.d0) & 
            +AGauss2*dexp(-((log(1.0d0+z)-zGauss2)/wGauss2)**2.d0))
    end if

    !!c	add the HeI part, using same T_0 and T_1 values
    Rdown_trip = a_trip/(sq_0*(1.d0+sq_0)**(1.0-b_trip))
    Rdown_trip = Rdown_trip/((1.d0+sq_1)**(1.d0+b_trip))
    Rup_trip = Rdown_trip*dexp(-h_P*C*L_He2St_ion/(k_B*Tmat))
    Rup_trip = Rup_trip*((CR*Tmat)**(1.5d0))*(4.d0/3.d0)
    !!c	last factor here is the statistical weight

    !!c       try to avoid "NaN" when x_He gets too small
    if ((x_He.lt.5.d-9) .or. (x_He.gt.0.980)) then 
       Heflag = 0
    else
       Heflag = Heswitch
    end if
    if (Heflag.eq.0)then		!use Peebles coeff. for He
       K_He = CK_He/Hz
    else	!for Heflag>0 		!use Sobolev escape probability
       tauHe_s = A2P_s*CK_He*3.d0*n_He*(1.d0-x_He)/Hz
       pHe_s = (1.d0 - dexp(-tauHe_s))/tauHe_s
       K_He = 1.d0/(A2P_s*pHe_s*3.d0*n_He*(1.d0-x_He))
       !!c	smoother criterion here from Antony Lewis & Chad Fendt
       if (((Heflag.eq.2).or.(Heflag.ge.5)).and.(x_H.lt.0.9999999d0))then
          !!c	use fitting formula for continuum opacity of H
          !!c	first get the Doppler width parameter
          Doppler = 2.D0*k_B*Tmat/(m_H*not4*C*C)
          Doppler = C*L_He_2p*dsqrt(Doppler)
          gamma_2Ps = 3.d0*A2P_s*ARGS_FHE*(1.d0-x_He)*C*C & 
               /(dsqrt(Pi)*sigma_He_2Ps*8.d0*Pi*Doppler*(1.d0-x_H)) & 
               /((C*L_He_2p)**2.d0)
          pb = 0.36d0  !value from KIV (2007)
          qb = b_He

          !!c	calculate AHcon, the value of A*p_(con,H) for H continuum opacity
          AHcon = A2P_s/(1.d0+pb*(gamma_2Ps**qb))
          K_He=1.d0/((A2P_s*pHe_s+AHcon)*3.d0*n_He*(1.d0-x_He))
       end if
       if (Heflag.ge.3) then		!include triplet effects
          tauHe_t = A2P_t*n_He*(1.d0-x_He)*3.d0
          tauHe_t = tauHe_t /(8.d0*Pi*Hz*L_He_2Pt**(3.d0))
          pHe_t = (1.d0 - dexp(-tauHe_t))/tauHe_t
          CL_PSt = h_P*C*(L_He_2Pt - L_He_2st)/k_B
	  if ((Heflag.eq.3) .or. (Heflag.eq.5).or.(x_H.gt.0.99999d0)) then
      !!c	no H cont. effect
             CfHe_t = A2P_t*pHe_t*dexp(-CL_PSt/Tmat)
             CfHe_t = CfHe_t/(Rup_trip+CfHe_t)	!"C" factor for triplets
	  else					!include H cont. effect
             Doppler = 2.d0*k_B*Tmat/(m_H*not4*C*C)
             Doppler = C*L_He_2Pt*dsqrt(Doppler)
             gamma_2Pt = 3.d0*A2P_t*ARGS_FHE*(1.d0-x_He)*C*C & 
                  /(dsqrt(Pi)*sigma_He_2Pt*8.d0*Pi*Doppler*(1.d0-x_H)) & 
                  /((C*L_He_2Pt)**2.d0)
             !!c	use the fitting parameters from KIV (2007) in this case
             pb = 0.66d0
             qb = 0.9d0
             AHcon = A2P_t/(1.d0+pb*gamma_2Pt**qb)/3.d0
             CfHe_t = (A2P_t*pHe_t+AHcon)*dexp(-CL_PSt/Tmat)
             CfHe_t = CfHe_t/(Rup_trip+CfHe_t)	!"C" factor for triplets
	  end if
       end if
    end if

    !!c	Estimates of Thomson scattering time and Hubble time
    timeTh=(1.d0/(CT*Trad**4))*(1.d0+x+ARGS_FHE)/x	!Thomson time
    timeH=2.d0/(3.d0*ARGS_H0*(1.d0+z)**1.5)		!Hubble time

    !!c	calculate the derivatives
    !!c	turn on H only for x_H<0.99, and use Saha derivative for 0.98<x_H<0.99
    !!c	(clunky, but seems to work)
    if (x_H.gt.0.99d0) then			!don't change at all
       f(1) = 0.d0
       !!cc	else if ((x_H.gt.0.98d0).and.(Heflag.eq.0)) then	!don't modify
    else if (x_H.gt.0.985d0) then		!use Saha rate for Hydrogen
       f(1) = (x*x_H*n*Rdown - Rup*(1.d0-x_H)*dexp(-CL/Tmat)) & 
            /(Hz*(1.d0+z))
       !!c	for interest, calculate the correction factor compared to Saha
       !!c	(without the fudge)
       factor=(1.d0 + K*Lambda*n*(1.d0-x_H)) & 
            /(Hz*(1.d0+z)*(1.d0+K*Lambda*n*(1.d0-x) & 
            +K*Rup*n*(1.d0-x)))
    else					!use full rate for H
       f(1) = ((x*x_H*n*Rdown - Rup*(1.d0-x_H)*dexp(-CL/Tmat)) & 
            *(1.d0 + K*Lambda*n*(1.d0-x_H))) & 
            /(Hz*(1.d0+z)*(1.d0/ARGS_FUDGE+K*Lambda*n*(1.d0-x_H)/ARGS_FUDGE & 
            +K*Rup*n*(1.d0-x_H)))
    end if
    !!c	turn off the He once it is small
    if (x_He.lt.1.d-15) then
       f(2)=0.d0
    else
       f(2) = ((x*x_He*n*Rdown_He  & 
            - Rup_He*(1.d0-x_He)*dexp(-CL_He/Tmat)) & 
            *(1.d0+ K_He*Lambda_He*n_He*(1.d0-x_He)*He_Boltz)) & 
            /(Hz*(1.d0+z) & 
            * (1.d0 + K_He*(Lambda_He+Rup_He)*n_He*(1.d0-x_He)*He_Boltz))
       !!c	Modification to HeI recombination including channel via triplets
       if (Heflag.ge.3) then
          f(2) = f(2)+ (x*x_He*n*Rdown_trip & 
               - (1.d0-x_He)*3.d0*Rup_trip*dexp(-h_P*C*L_He_2st/(k_B*Tmat))) & 
               *CfHe_t/(Hz*(1.d0+z))
       end if
    end if

    !!c	follow the matter temperature once it has a chance of diverging

    if (timeTh.lt.H_frac*timeH) then
       !!c		f(3)=Tmat/(1.d0+z)	!Tmat follows Trad
       !!c	additional term to smooth transition to Tmat evolution,
       !!c	(suggested by Adam Moss)
       epsilon = Hz*(1.d0+x+ARGS_FHE)/(CT*Trad**3*x)
       f(3) = bg%Tcmb() & 
            + epsilon*((1.d0+ARGS_FHE)/(1.d0+ARGS_FHE+x))*((f(1)+ARGS_FHE*f(2))/x) & 
            - epsilon* dHdz/Hz + 3.0d0*epsilon/(1.d0+z) 
    else
       f(3)= CT * (Trad**4) * x / (1.d0+x+ARGS_FHE) & 
            * (Tmat-Trad) / (Hz*(1.d0+z)) + 2.d0*Tmat/(1.d0+z)
    end if

    return

  end subroutine coop_recfast_ion


end module coop_recfast_mod
