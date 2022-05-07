#define USER_PARAM(i)  args%r(i)
!! you are allowed to use up to 10 parameters for primordial power 

!!Input:
!!kbykpiv =  k / k_{pivot}
!!cosmology
!!args (contains the user defined parameters)
!!Output:
!!ps = dimensionless scalar  primordial  power spectrum
!!pt = dimensionless tensor primordial power spectrum

!!100 is a sample code (not used anywhere); do not change it; the user-defined subroutines are from 101 to 110.
subroutine coop_user_defined_primordial_power_100(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
  COOP_REAL::nt
#define SCALAR_AS USER_PARAM(1)
#define SCALAR_NS USER_PARAM(2)
#define TENSOR_R  USER_PARAM(3)
  if(cosmology%inflation_consistency)then
     nt = - TENSOR_R/8.d0
  else
     nt = 0.d0
  endif
  ps =  SCALAR_AS * kbykpiv ** (SCALAR_NS - 1.d0)
  pt =  (SCALAR_AS * TENSOR_R) * kbykpiv ** nt
#undef SCALAR_AS
#undef SCALAR_NS
#undef TENSOR_R
end subroutine coop_user_defined_primordial_power_100


!!!!!!!!!!!!!!!!!!!!!! sin (sharp feature) 
subroutine coop_user_defined_primordial_power_101(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
!!write your code here
#define AS USER_PARAM(1)
#define NS USER_PARAM(2)
#define CF USER_PARAM(3)
#define KF USER_PARAM(4)
#define PHI USER_PARAM(5)
ps = AS * kbykpiv ** (NS - 1.d0) * (1.d0 + CF * sin(2 * kbykpiv * 0.05 / KF + PHI))
pt = 0.d0
!!to use this model; put "pp_genre = 101" in your ini file
!!stop "primordial power 101 has not been defined."
#undef AS
#undef NS
#undef CF
#undef KF
#undef PHI
end subroutine coop_user_defined_primordial_power_101


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! sin-log (resonance signal) 
subroutine coop_user_defined_primordial_power_102(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
!!write your code here
#define AS USER_PARAM(1)
#define NS USER_PARAM(2)
#define CF USER_PARAM(3)
#define OMEGA USER_PARAM(4)
#define PHI USER_PARAM(5)
ps = AS * kbykpiv ** (NS - 1.d0) * (1.d0 + CF * sin(OMEGA * log(2 * kbykpiv * 0.05) + PHI))
pt = 0.d0
!!to use this model; put "pp_genre = 102" in your ini file
!!stop "primordial power 102 has not been defined."
#undef AS
#undef NS
#undef CF
#undef OMEGA
#undef PHI
end subroutine coop_user_defined_primordial_power_102


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! full standard clock signal
subroutine coop_user_defined_primordial_power_103(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  real:: clock,delta, ksharp, ka, kb, sharp,ALPHA,BETA, OMEGA
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
!!write your code here
#define AS USER_PARAM(1)
#define NS USER_PARAM(2)
#define CF USER_PARAM(3)
#define KR USER_PARAM(4)
          ksharp=KR/(1.05*OMEGA)
	ALPHA=0.75*3.14159
	BETA=0.55*3.14159
	OMEGA=30
 clock = CF*(2.0*0.05*kbykpiv/KR)**(-3.0/2.0)*(sin(OMEGA*log(2.0*0.05*kbykpiv/KR)+ALPHA))
	  sharp = CF*(7.0*10.0**(-4.0)*(2.0*0.05*kbykpiv/ksharp)**2.0+0.5)*(cos(2.0*0.05*kbykpiv/ksharp+BETA))
	  ka=(67.0/140.0)*KR
	  kb=(96.0/140.0)*KR
	if (kbykpiv<0) then
	delta=0
	else if (kbykpiv<ka/0.05) then
	delta=sharp
	else if (kbykpiv>=ka/0.05 .and. kbykpiv<kb/0.05) then
	delta=(70.0/65.0)*clock
	else if (kbykpiv>=kb/0.05) then
	delta=(95.0/65.0)*clock
	end if
ps = AS * kbykpiv ** (NS - 1.d0) * (1.d0 + delta)
pt = 0.d0
!!to use this model; put "pp_genre = 103" in your ini file
!!stop "primordial power 103 has not been defined."
#undef AS
#undef NS
#undef CF
#undef ALPHA
#undef BETA
#undef OMEGA
#undef KR
#undef ka
#undef kb
#undef sharp
#undef ksharp
#undef delta
#undef clock
end subroutine coop_user_defined_primordial_power_103



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! clock signal for general scenario (arbitrary p)
subroutine coop_user_defined_primordial_power_104(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
  real::  clock_gen, delta
!!write your code here
#define AS USER_PARAM(1)
#define NS USER_PARAM(2)
#define CF USER_PARAM(3)
#define kr USER_PARAM(4)
#define pp USER_PARAM(5)
#define OMEGA USER_PARAM(6)
#define PHI USER_PARAM(7)
          
         clock_gen=CF*(2.0*kbykpiv*0.05/kr)**(-3.0/2.0+1/(2.0*pp))*(sin(pp*OMEGA/2*(2.0*kbykpiv*0.05/kr)**(1/pp)+PHI))
	if (pp<1 .and. pp>0) then
		if (kbykpiv<kr/(2*0.05) .and. kbykpiv > kr/(OMEGA*0.05)) then
			delta=clock_gen
		else
			delta=0
		end if
	else 
		if (kbykpiv>kr/(2*0.05)) then
			delta=clock_gen
		else
			delta=0
		end if
	end if
ps = AS * kbykpiv ** (NS - 1.d0) * (1.d0 + delta)
pt = 0.d0
!!to use this model; put "pp_genre = 104" in your ini file
!!stop "primordial power 104 has not been defined."
#undef AS
#undef NS
#undef CF
#undef kr
#undef pp
#undef OMEGA
#undef PHI
#undef clock_gen
#undef delta
#undef clock
end subroutine coop_user_defined_primordial_power_104



!! sharp feature with Cora's envelop
subroutine coop_user_defined_primordial_power_105(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
  real:: CalD, xx
!!write your code here
#define AS USER_PARAM(1)
#define NS USER_PARAM(2)
#define CF USER_PARAM(3)
#define etaF USER_PARAM(4)
#define xd USER_PARAM(5)
xx = 0.05 * kbykpiv * etaF
CalD = xx/(xd* sinh(xx/xd)) * ((-3+9/xx**2) * cos( 2* xx) + (15 - 9/ xx**2) *sin(2* xx)/(2* xx))
ps = AS * kbykpiv ** (NS - 1.d0) * exp(CF* CalD)
pt = 0.d0
!!to use this model; put "pp_genre = 105" in your ini file
!!stop "primordial power 105 has not been defined."
#undef AS
#undef NS
#undef CF
#undef etaF
#undef xx
#undef CalD
#undef xd
end subroutine coop_user_defined_primordial_power_105


!! localized feature: particle production
subroutine coop_user_defined_primordial_power_106(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
  real:: Pi
!!write your code here
#define AS USER_PARAM(1)
#define NS USER_PARAM(2)
#define CF USER_PARAM(3)
#define KF USER_PARAM(4)
Pi = 3.1415927
ps =  AS * kbykpiv ** (NS - 1.d0) * (1+ CF* (kbykpiv*0.05/KF)**3 * exp( - Pi/2 * (kbykpiv*0.05/KF)**2 ))
pt = 0.d0
!!to use this model; put "pp_genre = 106" in your ini file
!!stop "primordial power 106 has not been defined."
#undef AS
#undef NS
#undef CF
#undef KF
#undef Pi
end subroutine coop_user_defined_primordial_power_106


!!!!!!!!!!!!!!!!!!!! SC with no clock signal!!!!!!!!!!!!!!!!!1
subroutine coop_user_defined_primordial_power_107(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  real:: delta, ksharp, ka, sharp,BETA, OMEGA
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
!!write your code here
#define AS USER_PARAM(1)
#define NS USER_PARAM(2)
#define CF USER_PARAM(3)
#define KR USER_PARAM(4)
          ksharp=KR/(1.05*OMEGA)
	BETA=0.55*3.14159
	OMEGA=30
	  sharp = CF*(7.0*10.0**(-4.0)*(2.0*0.05*kbykpiv/ksharp)**2.0+0.5)*(cos(2.0*0.05*kbykpiv/ksharp+BETA))
	  ka=(67.0/140.0)*KR
	if (kbykpiv<0) then
	delta=0
	else if (kbykpiv<ka/0.05) then
	delta=sharp
	else if (kbykpiv>=ka/0.05) then
	delta=0
	end if
ps = AS * kbykpiv ** (NS - 1.d0) * (1.d0 + delta)
pt = 0.d0
!!to use this model; put "pp_genre = 107" in your ini file
!!stop "primordial power 103 has not been defined."
#undef AS
#undef NS
#undef CF
#undef BETA
#undef OMEGA
#undef KR
#undef ka
#undef sharp
#undef ksharp
#undef delta
end subroutine coop_user_defined_primordial_power_107


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! clock signal for general scenario (arbitrary p). 
!!!! This is the same as clock signal in _104 but we defined here: C_eff = C*(2 k_pivot/k_r)^{-3/2+1/(2p)} and OMEGA_eff = Omega (2 k_pivot /k_r)^{1/p};
subroutine coop_user_defined_primordial_power_108(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
  real::  clock_gen, delta
!!write your code here
#define AS USER_PARAM(1)
#define NS USER_PARAM(2)
#define CF USER_PARAM(3)
#define kr USER_PARAM(4)
#define pp USER_PARAM(5)
#define OMEGA USER_PARAM(6)
#define PHI USER_PARAM(7)
          
  clock_gen=CF*(kbykpiv)**(-3.0/2.0+1/(2.0*pp))* (sin(pp*OMEGA*(kbykpiv)**(1/pp)/2+PHI)) 
  !         clock_gen=CF*(2.0*kbykpiv*0.05/kr)**(-3.0/2.0+1/(2.0*pp))* (sin(pp*OMEGA/2*(2.0*kbykpiv*0.05/kr)**(1/pp)+PHI))
  if (pp<1 .and. pp>0) then

     if (kbykpiv<kr/(2*0.05) .and. kbykpiv > kr*(2*0.05/kr)**(1/pp)/(OMEGA*0.05)) then
        delta= clock_gen
     else
        delta=0
     end if
  else 
     if (kbykpiv>kr/(2*0.05)) then
        delta=clock_gen
     else
        delta=0
     end if
  end if
ps = AS * kbykpiv ** (NS - 1.d0) * (1.d0 + delta)
pt = 0.d0
!!to use this model; put "pp_genre = 108" in your ini file
!!stop "primordial power 104 has not been defined."
#undef AS
#undef NS
#undef CF
#undef kr
#undef pp
#undef OMEGA
#undef PHI
#undef clock_gen
#undef delta
#undef clock
end subroutine coop_user_defined_primordial_power_108


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! clock signal for general scenario (arbitrary p). 
!!!! This is the same as clock signal in _104 but we defined here: C_eff = C*(2 k_pivot/k_r)^{-3/2+1/(2p)} and OMEGA_eff = Omega (2 k_pivot /k_r)^{1/p} and phi_eff=phi+p Omega_eff/2. This is appropriate for the inflationary clock signal
subroutine coop_user_defined_primordial_power_109(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
  real::  clock_gen, delta
!!write your code here
#define AS USER_PARAM(1)
#define NS USER_PARAM(2)
#define CF USER_PARAM(3)
#define kr USER_PARAM(4)
#define pp USER_PARAM(5)
#define OMEGA USER_PARAM(6)
#define PHI USER_PARAM(7)
          
  	clock_gen=CF*(kbykpiv)**(-3.0/2.0+1/(2.0*pp))* (sin(pp*OMEGA*(kbykpiv)**(1/pp)/2+PHI-pp*OMEGA/2)) 
	!clock_gen=CF*(kbykpiv)**(-3.0/2.0+1/(2.0*pp))* (sin(pp*OMEGA*(kbykpiv)**(1/pp)/2+PHI-pp*OMEGA/2.0 + OMEGA*log(0.05)/2.0 )) 
  if (pp<1 .and. pp>0) then

     if (kbykpiv<kr/(2*0.05) .and. kbykpiv > kr*(2*0.05/kr)**(1/pp)/(OMEGA*0.05)) then
        delta= clock_gen
     else
        delta=0
     end if
  else 
     if (kbykpiv>kr/(2*0.05)) then
        delta=clock_gen
     else
        delta=0
     end if
  end if
ps = AS * kbykpiv ** (NS - 1.d0) * (1.d0 + delta)
pt = 0.d0
!!to use this model; put "pp_genre = 109" in your ini file
!!stop "primordial power 109 has not been defined."
#undef AS
#undef NS
#undef CF
#undef kr
#undef pp
#undef OMEGA
#undef PHI
#undef clock_gen
#undef delta
#undef clock
end subroutine coop_user_defined_primordial_power_109


!!!! clock signal for infaltionary scenario (resonance signal with (k^{-3/2}) envelop)
subroutine coop_user_defined_primordial_power_110(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
  real::  clock_gen, delta
!!write your code here
#define AS USER_PARAM(1)
#define NS USER_PARAM(2)
#define CF USER_PARAM(3)
#define kr USER_PARAM(4)
#define OMEGA USER_PARAM(5)
#define PHI USER_PARAM(6)
          
  	clock_gen=CF*(2*kbykpiv*0.05/kr)**(-3.0/2.0)* (sin(OMEGA*log(2*kbykpiv*0.05)/2+PHI)) 

     if (kbykpiv>kr/(2*0.05)) then
        delta=clock_gen
     else
        delta=0
     end if
ps = AS * kbykpiv ** (NS - 1.d0) * (1.d0 + delta)
pt = 0.d0
!!to use this model; put "pp_genre = 110" in your ini file
!!stop "primordial power 110 has not been defined."
#undef AS
#undef NS
#undef CF
#undef kr
#undef OMEGA
#undef PHI
#undef clock_gen
#undef delta
#undef clock
end subroutine coop_user_defined_primordial_power_110





#undef USER_PARAM
