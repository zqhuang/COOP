  if(this%mask(i, inow) .eq. 1)then
     lapln =  dot_product(this%lnphi(i-1:i+1, inow), this%lapc(:, i)) + (this%lnphi(i+1, inow)-this%lnphi(i-1, inow))**2/fourdr2
     phi = exp(this%lnphi(i, inow))
     err  = phi*lapln - ((exp(this%lnrho(i, inow))-this%rho0)/this%a(inow)**3 - this%deltaRofphi(exp(this%lnphi(i, inow)), inow)) * a2by3
     this%lnphi(i, inow) = this%lnphi(i, inow) + min(max( err/((-this%lapc(0, i) - lapln + this%m2ofphi(phi) * a2by3)*phi), -0.05d0), 0.05d0)
     s = s + abs(err)
  endif
