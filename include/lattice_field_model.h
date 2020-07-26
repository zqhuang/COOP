!!in this file you define the potential and its derivatives
#define PHI fields(1)
#define CHI fields(2)  

  function coop_lattice_fields_V(fields) result(V)
    COOP_REAL::fields(:), V
    COOP_REAL,parameter::mphi = 1.191e-5*coop_lattice_Mp
    V = (0.75d0*(mphi*coop_lattice_Mp)**2)*(1.d0-exp(-sqrt(2.d0/3.d0)*PHI/coop_lattice_Mp))**2 * (1.d0 + alpha*exp(-(PHI-phipk)**2/2.d0/mu**2))
    !(mphi*PHI)**2/2.d0*(1.d0 + alpha*exp(-(PHI-phipk)**2/2.d0/mu**2))
  end function coop_lattice_fields_V

!!first derivatives  
  function coop_lattice_fields_dVdphi(fields) result(dVdphi)
    COOP_REAL,parameter::step=coop_lattice_Mp/(2**21)
    COOP_REAL::fields(:), dVdphi(size(fields)), tmp(size(fields)), Vup, Vdown
    COOP_INT::i
    tmp = fields
    do i=1, size(fields)
       tmp(i) = fields(i) + step
       Vup = coop_lattice_fields_V(tmp)
       tmp(i) = fields(i) - step
       Vdown = coop_lattice_fields_V(tmp)
       dVdphi(i) = (Vup - Vdown)/(2.d0*step)       
       tmp(i) = fields(i)
    enddo
  end function coop_lattice_fields_dVdphi

!!second derivatives   
  function coop_lattice_fields_d2Vdphi2(fields) result(d2Vdphi2)
    COOP_REAL,parameter::step=coop_lattice_Mp/(2**20)
    COOP_REAL::fields(:)
    COOP_REAL::d2Vdphi2(size(fields), size(fields)), tmp(size(fields)), Vpp, Vdd, Vpd, Vdp, Vc
    COOP_INT::i,j
    tmp = fields
    Vc = coop_lattice_fields_V(fields)    
    do i=1, size(fields)
       tmp(i) = fields(i) + step
       Vpp = coop_lattice_fields_V(tmp)
       tmp(i) = fields(i) - step
       Vdd = coop_lattice_fields_V(tmp)
       d2Vdphi2(i,i) = (Vpp+Vdd-2.d0*Vc)/(step**2)
       tmp(i) = fields(i)
       do j=1, i-1
          tmp(i) = fields(i) + step
          tmp(j) = fields(j) + step
          Vpp = coop_lattice_fields_V(tmp)
          tmp(i) = fields(i) - step
          tmp(j) = fields(j) - step
          Vdd = coop_lattice_fields_V(tmp)
          tmp(i) = fields(i) + step
          tmp(j) = fields(j) - step
          Vpd = coop_lattice_fields_V(tmp)
          tmp(i) = fields(i) - step
          tmp(j) = fields(j) + step
          Vdp = coop_lattice_fields_V(tmp)
          d2Vdphi2(i, j) = (Vpp + Vdd - Vpd - Vdp)/(2.d0*step**2)
          d2Vdphi2(j, i) = d2Vdphi2(i, j)
          tmp(i) = fields(i)
          tmp(j) = fields(j)
       enddo
    enddo

  end function coop_lattice_fields_d2Vdphi2
  
#undef PHI
#undef CHI  
