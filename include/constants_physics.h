!!physical constants
COOP_REAL, parameter::coop_SI_c = 299792458.d0
COOP_REAL, parameter::coop_SI_kB =  1.3806504d-23
COOP_REAL, parameter::coop_SI_h = 6.626068d-34
COOP_REAL, parameter::coop_SI_hbar = coop_SI_h / (2.d0*coop_pi)
COOP_REAL, parameter::coop_SI_G = 6.67428d-11
COOP_REAL, parameter::coop_SI_alpha = 1.d0/137.035d0
!!units
!!angle
COOP_REAL, parameter::coop_SI_degree = coop_pi/180.d0
COOP_REAL, parameter::coop_SI_arcmin = coop_SI_degree/60.d0
COOP_REAL, parameter::coop_SI_arcsec = coop_SI_arcmin / 60.d0
COOP_REAL, parameter::coop_fullsky_degrees = 4.d0*coop_pi / coop_SI_degree**2
COOP_REAL, parameter::coop_sigma_by_fwhm = 1.d0/sqrt(8.d0*coop_ln2)
!!time
doubleprecision, parameter::coop_SI_PlanckTime = 5.39124d-44
COOP_REAL, parameter::coop_SI_day = 3600.d0*24.d0
COOP_REAL, parameter::coop_SI_Year =  coop_SI_day*365.2422d0
COOP_REAL, parameter::coop_SI_kyr =  coop_SI_Year * 1.d3
COOP_REAL, parameter::coop_SI_Myr =  coop_SI_Year * 1.d6
COOP_REAL, parameter::coop_SI_Gyr =  coop_SI_Year * 1.d9
!!length
doubleprecision, parameter::coop_SI_PlanckLength = 1.616252d-35
COOP_REAL, parameter::coop_SI_pc = 3.08568025d16
COOP_REAL, parameter::coop_SI_kpc = coop_SI_pc * 1.d3
COOP_REAL, parameter::coop_SI_Mpc = coop_SI_pc * 1.d6
COOP_REAL, parameter::coop_SI_Gpc = coop_SI_pc * 1.d9
COOP_REAL, parameter::coop_SI_lyr = coop_SI_c * coop_SI_year
!!energy
COOP_REAL, parameter::coop_SI_eV = 1.60217648740d-19
COOP_REAL, parameter::coop_SI_keV = coop_SI_eV * 1.d3
COOP_REAL, parameter::coop_SI_MeV = coop_SI_eV * 1.d6
COOP_REAL, parameter::coop_SI_GeV = coop_SI_eV * 1.d9
doubleprecision, parameter::coop_SI_PlanckEnergy = 1.2209d19 * coop_SI_GeV
COOP_REAL, parameter::coop_SI_H_bind = 13.605698 * coop_SI_eV
!!mass
doubleprecision, parameter::coop_SI_PlanckMass = coop_SI_PlanckEnergy/coop_SI_c**2
doubleprecision, parameter::coop_SI_Reduced_PlanckMass = coop_SI_PlanckMass/sqrt(8.d0*coop_pi)
COOP_REAL, parameter::coop_SI_Msun =  1.98892d30
COOP_REAL, parameter::coop_SI_Mearth = 5.97219d24
COOP_REAL, parameter::coop_SI_atomic_mass_unit =  1.66053878283d-27
COOP_REAL, parameter::coop_SI_m_e = coop_SI_atomic_mass_unit / 1822.8884845d0
COOP_REAL, parameter::coop_SI_m_p = 1.00727646681290d0*coop_SI_atomic_mass_unit
COOP_REAL, parameter::coop_SI_m_n = 1.0086649160043*coop_SI_atomic_mass_unit
COOP_REAL, parameter::coop_SI_m_H =  coop_SI_m_e + coop_SI_m_p - coop_SI_H_bind/coop_SI_c**2
COOP_REAL, parameter::coop_m_He_by_m_H = 3.9715d0
COOP_REAL, parameter::coop_SI_barssc0 = coop_SI_kB /coop_SI_m_p/coop_SI_c**2
!!temperature
doubleprecision, parameter::coop_SI_PlanckTemperature = coop_SI_PlanckEnergy / coop_SI_kB
!!blackbody
COOP_REAL, parameter::coop_SI_blackbody_alpha = 3.7828737060437682d-16  !!blackbody radiation density = alpha * (total g) * (Temperature in Kelvin) ^ 4  [g=1 for each spin degree of boson, g=7/8 for each spin degree of fermion; for photon g_total = 2; for neutrinos g_total = 7/8 * number of species] 
COOP_REAL, parameter::coop_SI_Stefan_Boltzmann = 5.6704d-8 !!radiation power per surface area pur Kelvin (W/m^2/K^4)
COOP_REAL, parameter::coop_SI_sigma_thomson = 6.6524616d-29 !!Thomson scattering cross section
!!cosmology
COOP_REAL, parameter::coop_SI_rhocritbyh2 = 1.8783467654345834d-26*coop_SI_c**2
COOP_REAL, parameter::coop_neutrinos_temperature_correction = (3.046d0/3.d0)**(1.d0/4.d0)
COOP_REAL, parameter::coop_SI_hbyH0 = coop_SI_Mpc /1.d5


