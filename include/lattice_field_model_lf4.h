!!in this file you define the potential and its derivatives
#define PHI phi(1)
#define CHI phi(2)  

  function coop_lattice_fields_V(phi) result(V)
    COOP_REAL::phi(:), V
    V =  lambda * ( (1.d0/4.d0) * PHI **2 + (g2byl/2.d0) * CHI **2 ) * PHI**2
  end function coop_lattice_fields_V

!!first derivatives  
  function coop_lattice_fields_dVdphi(phi) result(dVdphi)
    COOP_REAL::phi(:), dVdphi(size(phi))
    dVdphi(1) = lambda*(PHI**2 + g2byl * CHI**2) * PHI
    dVdphi(2) = (lambda*g2byl) * CHI * PHI**2
  end function coop_lattice_fields_dVdphi

!!second derivatives   
  function coop_lattice_fields_d2Vdphi2(phi) result(d2Vdphi2)
    COOP_REAL::phi(:)
    COOP_REAL::d2Vdphi2(size(phi), size(phi))
    d2Vdphi2(1,1) = (3.d0*PHI**2 + g2byl * CHI**2) * lambda
    d2Vdphi2(2,2) = (g2byl*lambda)*PHI**2
    d2Vdphi2(1,2) = (2.d0*g2byl*lambda)*PHI*CHI
    d2Vdphi2(2,1) = d2Vdphi2(1,2)
  end function coop_lattice_fields_d2Vdphi2
  
#undef PHI
#undef CHI  
