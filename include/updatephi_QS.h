!!$             this%lnphi(i, inow) = this%lnphi(i, ill) + this%lnphip(i, il)*this%twodtau
!!$             this%lnphip(i, inow) = this%lnphip(i, ill) + ( &
!!$                  (dot_product(this%lnphi(i-1:i+1, il), this%lapc(:, i)) + (this%lnphi(i+1, il)-this%lnphi(i-1, il))**2/fourdr2) &
!!$                  + (this%deltaRofphi(exp(this%lnphi(i, il)), il) - (exp(this%lnrho(i, il))-this%rho0)/this%a(il)**3 ) * this%a(il)**2/3.d0 *exp(-this%lnphi(i, il))  &
!!$                  - 2.d0*this%H(il)*this%lnphip(i, il)-this%lnphip(i, il)**2)*this%twodtau

  if(this%mask(i, inow) .eq. 1)then
     lapln =  dot_product(this%lnphi(i-1:i+1, inow), this%lapc(:, i)) + (this%lnphi(i+1, inow)-this%lnphi(i-1, inow))**2/fourdr2
     phi = exp(this%lnphi(i, inow))
     err  = phi*lapln - ((exp(this%lnrho(i, inow))-this%rho0)/this%a(inow)**3 - this%deltaRofphi(phi, inow)) * a2by3
     this%lnphi(i, inow) = this%lnphi(i, inow) + min(max( err/((-this%lapc(0, i) - lapln + this%m2ofphi(phi) * a2by3)*phi), -0.05d0), 0.05d0)
     s = s + abs(err)
  endif
