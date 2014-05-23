module coop_constants
  implicit none

#include "constants.h"
  integer,parameter::coop_real_length = 8
  integer,parameter::coop_short_string_length = 32
  integer,parameter::coop_string_length = 256
  integer,parameter::coop_long_string_length = 8192

  Character,Parameter::coop_backslash = Char(92)
  character,parameter::coop_slash = Char(47)
  Character,Parameter::coop_backspace = Char(8)
  Character,Parameter::coop_tab = Char(9)
  Character,Parameter::coop_newline = Char(10)
  Character,Parameter::coop_vertical_tab = Char(11)
  Character,Parameter::coop_newpage = Char(12)
  Character,Parameter::coop_carriage_return = Char(13)
  real(dl),parameter:: coop_pi = 3.14159265358979323846264338327950288d0
  real(dl),parameter:: coop_pio2 = coop_pi/2._dl
  real(dl),parameter:: coop_pio4 = coop_pi/4._dl
  real(dl),parameter:: coop_ln2 = 0.6931471805599453094172321d0 
  real(dl),parameter:: coop_ln10 = 2.302585092994045684017991d0
  real(dl),PARAMETER:: coop_LnPi=1.144729885849400174143427351d0
  real(dl),PARAMETER:: coop_LogPi=0.4971498726941338543512682883d0
  real(dl),parameter:: coop_sqrt2 = 1.4142135623730950488016887d0
  real(dl),PARAMETER:: coop_sqrt3 = 1.73205080756887729352744634d0
  real(dl),parameter:: coop_sqrt5 = 2.236067977499789696409174d0
  real(dl),parameter:: coop_sqrt6 = coop_sqrt2*coop_sqrt3
  real(dl),parameter:: coop_sqrt7 = 2.645751311064590590502d0
  real(dl),parameter:: coop_sqrt8 = 2.d0*coop_sqrt2
  real(dl),parameter:: coop_sqrt10 = coop_sqrt2 * coop_sqrt5
  real(dl),parameter:: coop_sqrt11 = 3.316624790355399849115d0
  real(dl),parameter:: coop_sqrt12 = 2.d0*coop_sqrt3
  real(dl),parameter:: coop_sqrt13 = 3.605551275463989293119d0
  real(dl),parameter:: coop_sqrt14 = coop_sqrt2*coop_sqrt7
  real(dl),parameter:: coop_sqrt15 = coop_sqrt3*coop_sqrt5
  real(dl),parameter:: coop_8pi = coop_pi*8._dl
  real(dl),parameter:: coop_4pi = coop_pi*4._dl
  real(dl),parameter:: coop_4piby3 = coop_4pi/3._dl
  real(dl),parameter:: coop_3by4pi = 3._dl/coop_4pi
  real(dl),parameter:: coop_2pi = coop_pi*2._dl
  real(dl),parameter:: coop_pi2 = coop_pi ** 2
  real(dl),parameter:: coop_pi4 = coop_pi2 ** 2
  real(dl),parameter:: coop_2pi2 = coop_pi2 * 2._dl
  real(dl),parameter:: coop_8pi3 = coop_2pi ** 3
  real(dl),parameter:: coop_7pi4by120 = coop_pi4 * (7.d0/120.d0)
  real(dl),parameter:: coop_third = 1._dl/3._dl
  real(dl),parameter:: coop_two_thirds = 2._dl/3._dl
  real(dl),parameter:: coop_four_thirds = 4._dl/3._dl
  real(dl),parameter:: coop_sqrtpi = 1.7724538509055160272981674833411_dl
  real(dl),parameter:: coop_EulerC=0.57721566490153286060651209_dl
  real(dl),parameter:: coop_Riemannzeta3 = 1.2020569031595942853997_dl
  real(dl),parameter:: coop_Riemannzeta5  = 1.0369277551433699263313_dl
  real(dl),parameter:: coop_Riemannzeta7  = 1.0083492773819228268397_dl
  real(dl),parameter:: coop_Riemannzeta9 = 1.00200839282608221441785_dl
  real(dl),PARAMETER:: coop_fullsky_degrees = 41252.96125   !!4pi/degree^2
  real(dl),parameter:: coop_degree = coop_pi/180.d0
  real(dl),parameter:: coop_arcmin = coop_degree/60.d0
  real(dl),parameter:: coop_sigmabyfwhm = 1.d0/sqrt(8.d0*coop_ln2)
  real(dl),parameter:: coop_arcsec = coop_arcmin/60.d0
  real(dl),parameter:: coop_chbyMpcH0 = 2997.92458d0
  real(dl),parameter :: coop_NewtonG = 6.67428e-11_dl
  real(dl),parameter:: coop_planck = 6.626068d-34
  real(dl),parameter:: coop_hbar = coop_planck/coop_2pi
  real(dl),parameter:: coop_boltzmann = 1.3806504d-23
  real(dl),parameter:: coop_K2GHz = coop_boltzmann/(coop_planck*1.d9) !!
  real(dl),parameter:: coop_Mpc_SI = 3.08568025e22_dl
  real(dl),parameter:: coop_Msun_SI = 1.98892e30_dl
  real(dl),parameter:: coop_SpeedOfLight_SI = 299792458._dl
  real(dl),parameter:: coop_c = coop_SpeedOfLight_SI
  real(dl),parameter:: coop_PlanckMass_SI = 2.17644e-8_dl
  real(dl),parameter:: coop_PlanckLength_SI =1.616252e-35_dl
  real(dl),parameter:: coop_PlanckTime_SI = 5.39124e-44_dl
  real(dl),parameter:: coop_PlanckEnergy_SI = coop_PlanckMass_SI * coop_SpeedOfLight_SI ** 2
  real(dl),parameter:: coop_PlanckTemperature_SI = coop_PlanckEnergy_SI / coop_boltzmann
  real(dl),parameter:: coop_Yr_SI = 3600._dl*24._dl*365.2422_dl !!time
  real(dl),parameter:: coop_Gyr_SI = coop_Yr_SI * 1.e9_dl !!time
  real(dl),parameter:: coop_Lyr_SI = 9.4605284e15_dl !!length
  real(dl),parameter:: coop_eV_SI = 1.60217648740e-19_dl
  real(dl),parameter:: coop_eVmass_SI = coop_eV_SI / coop_SpeedOfLight_SI ** 2
  real(dl),parameter:: coop_eVlength_SI = coop_SpeedOfLight_SI/ (coop_eV_SI /coop_hbar)
  real(dl),parameter:: coop_MeV_SI = coop_eV_SI * 1.e6_dl
  real(dl),parameter:: coop_MeVmass_SI = coop_MeV_SI / coop_SpeedOfLight_SI ** 2
  real(dl),parameter:: coop_GeV_SI = coop_eV_SI * 1.e9_dl
  real(dl),parameter:: coop_GeVmass_SI = coop_GeV_SI / coop_SpeedOfLight_SI ** 2

  real(dl),parameter:: coop_fine_structure = 1._dl/137.035_dl
  real(dl),parameter:: coop_atomic_mass_unit_SI = 1.66053878283e-27_dl

  real(dl),parameter:: coop_electron_mass_SI = coop_atomic_mass_unit_SI / 1822.8884845_dl !!0.51099892811 * coop_MeVmass_SI

  real(dl),parameter:: coop_proton_mass_SI = 1.00727646681290 * coop_atomic_mass_unit_SI

  real(dl),parameter:: coop_E_hydrogen = 13.605698 * coop_eV_SI 

  real(dl),parameter:: coop_Hydrogen_mass_SI = coop_electron_mass_SI + coop_proton_mass_SI - coop_E_Hydrogen/coop_SpeedOfLight_SI**2

  real(dl),parameter:: coop_massratio_He_H =  3.9715_dl

  real(dl),parameter:: coop_rhocritbyh2_SI = 3._dl*(1.e5_dl/coop_SpeedOfLight_SI)**2*(coop_PlanckMass_SI/coop_Mpc_SI**2/coop_PlanckLength_SI)/coop_8pi !!1.878e-26 !!mass/volume, note that this is NOT energy per volume
  real(dl),parameter:: coop_hbyH0_SI = coop_Mpc_SI / 1.e5_dl !!Cosmic age unit
  real(dl),parameter:: coop_hbyH0_Gyr = coop_hbyH0_SI / coop_Gyr_SI
  real(dl),parameter:: coop_chbyH0_Mpc = coop_chbyMpcH0  !!alias
  real(dl),parameter:: coop_chbyH0_SI = coop_chbyH0_Mpc * coop_Mpc_SI  !!alias
  real(dl),parameter:: coop_H0byh_SI = 1./coop_hbyH0_SI
  real(dl),parameter:: coop_H0byh_eV = coop_hbar/coop_eV_SI * coop_H0byh_SI

  real(dl),parameter:: coop_blackbody_alpha_SI = (coop_pi2/30._dl)*(coop_boltzmann) / (coop_PlanckTemperature_SI * coop_PlanckLength_SI)**3  !!blackbody radiation density = alpha * (total g) * (Temperature in Kelvin) ^ 4  [g=1 for each spin degree of boson, g=7/8 for each spin degree of fermion; for photon g_total = 2; for neutrinos g_total = 7/8 * number of species] (result is in SI unit J/m^3)
  real(dl),parameter:: coop_Stefan_Boltzmann = 5.670400e-8_dl !!W/m^2/K^4

  real(dl),parameter:: coop_sigma_thomson_SI = 6.6524616e-29_dl

  real(dl), parameter :: coop_arad_SI = (coop_pi**5 * 8./15. ) * coop_boltzmann * (coop_boltzmann/(coop_SpeedOfLight_SI * coop_Planck)) **3  
    !7.565914e-16_dl !radiation constant for u=aT^4
  !! = 2.*coop_blackboday_alpha_SI (since photon has two spin dof)

  real(dl), parameter :: coop_ComptonCT = (8.d0/3.d0) * coop_sigma_thomson_SI/(coop_electron_mass_SI * coop_SpeedOfLight_SI**2) * coop_arad_SI  * coop_chbyH0_SI  !! ch/H_0 * 8/3 alpha * sigma_T / (m_e c^2)
  real(dl),parameter::coop_barssc0 = coop_boltzmann / coop_hydrogen_mass_SI / coop_SpeedOfLight_SI ** 2


  real(dl),parameter:: coop_rhocritbyh2_MsunbyMpc3 = 3._dl*(coop_PlanckMass_SI/coop_Msun_SI)*(coop_Mpc_SI/coop_PlanckLength_SI) * (1.e5_dl/coop_SpeedOfLight_SI)**2 / coop_8pi !!2.77467e11
  real(dl),parameter::coop_deltac_EdS = (3._dl/5._dl)*(3._dl/2._dl*coop_pi)**(2._dl/3._dl) !!spherical collapse critical density contrast


  

end module coop_constants
