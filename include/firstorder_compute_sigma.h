
  function coop_cosmology_firstorder_sigma_TopHat_R(this, z, R) result(sigma)
    class(coop_cosmology_firstorder)::this    
    COOP_INT, parameter::nk = 800
    COOP_INT ik
    COOP_REAL lnk(nk), k(nk), Pk(nk), dlnk
    COOP_REAL z, R, sigma
    call coop_set_uniform(nk, lnk, min(1.d0, -log(R)-2.d0), -log(R) + 3.25d0)
    k = exp(lnk)
    dlnk = lnk(2)-lnk(1)
    call this%get_matter_power(z, nk, k, pk)
    !$omp parallel do
    do ik = 1, nk
       pk(ik) = pk(ik)*FT_tophat(k(ik)*R)**2
    enddo
    !$omp end parallel do
    sigma = sqrt((sum(pk)*dlnk + (1.d0/3.d0-dlnk/2.d0)*pk(1)))*3.d0
  contains

    function FT_tophat(kR) 
      !!1/3 * Fourier transformation of a tophat function in a sphere with radius R
      COOP_REAL kR, FT_tophat
      if(kR .gt. 0.02d0)then
         FT_tophat = (dsin(kR)/kR -  dcos(kR))/ kR**2 
      else
         FT_tophat = (10.d0 - kR**2*(1.d0 - kR**2/28.d0))/30.d0
      endif
    end function FT_tophat
    
  end function coop_cosmology_firstorder_sigma_TopHat_R


  !!d1 = d ln (sigma^2) / d ln R
  !!d2 = d^2 ln (sigma^2) /d (ln R)^2
  subroutine coop_cosmology_firstorder_sigma_Gaussian_R_with_dervs(this, z, R, sigma, d1, d2)
    class(coop_cosmology_firstorder)::this    
    COOP_INT, parameter::nk = 1500
    COOP_INT ik
    COOP_REAL lnk(nk), k(nk), Pk(nk), Pk1(nk), Pk2(nk), dlnk, x2
    COOP_REAL z, R, sigma, d1, d2, sum1, sum2, sum3
    call coop_set_uniform(nk, lnk, min(0.5d0, -log(R)-1.5d0), -log(R)+1.5d0)
    k = exp(lnk)
    dlnk = lnk(2)-lnk(1)
    call this%get_matter_power(z, nk, k, pk)
    !$omp parallel do private(x2)
    do ik = 1, nk
       x2 = (k(ik)*R)**2
       pk(ik) = pk(ik)*exp(-x2)
       Pk1(ik) = pk(ik)*x2
       Pk2(ik) = pk(ik)*x2*(1.d0-x2)
    enddo
    !$omp end parallel do
    sum1 = (sum(pk)*dlnk + pk(1)*(1.d0/3.d0 - dlnk/2.d0))
    sum2 = (sum(pk1)*dlnk + pk1(1)*(1.d0/5.d0 - dlnk/2.d0))*2.d0
    sum3 = (sum(pk2)*dlnk + pk2(1)*(1.d0/5.d0 - dlnk/2.d0))*4.d0
    sigma = sqrt(sum1)
    d1 = -sum2/sum1
    d2 = -(sum2/sum1)**2 - sum3/sum1
  end subroutine coop_cosmology_firstorder_sigma_Gaussian_R_with_dervs


  function coop_cosmology_firstorder_sigma_Gaussian_R(this, z, R) result(sigma)
    class(coop_cosmology_firstorder)::this    
    COOP_INT, parameter::nk = 1800
    COOP_INT ik
    COOP_REAL lnk(nk), k(nk), Pk(nk), dlnk
    COOP_REAL z, R, sigma
    call coop_set_uniform(nk, lnk, min(0.5d0, -log(R)-1.5d0), -log(R)+1.5d0)
    k = exp(lnk)
    dlnk = lnk(2)-lnk(1)
    call this%get_matter_power(z, nk, k, pk)
    !$omp parallel do
    do ik = 1, nk
       pk(ik) = pk(ik)*exp(-(k(ik)*R)**2)
    enddo
    !$omp end parallel do
    sigma = sqrt(sum(pk)*dlnk + pk(1)*(1.d0/3.d0 - dlnk/2.d0))
  end function coop_cosmology_firstorder_sigma_Gaussian_R


  !!quick estimate
  function coop_cosmology_firstorder_sigma_Gaussian_R_quick(this, z, R) result(sigma)
    class(coop_cosmology_firstorder)::this    
    COOP_INT, parameter::nk = 600
    COOP_INT ik
    COOP_REAL lnk(nk), k(nk), Pk(nk), dlnk
    COOP_REAL z, R, sigma
    call coop_set_uniform(nk, lnk, min(1.5d0, -log(R)-1.d0), -log(R)+1.25d0)
    k = exp(lnk)
    dlnk = lnk(2)-lnk(1)
    call this%get_matter_power(z, nk, k, pk)
    !$omp parallel do
    do ik = 1, nk
       pk(ik) = pk(ik)*exp(-(k(ik)*R)**2)
    enddo
    !$omp end parallel do
    sigma = sqrt((sum(pk)*dlnk + pk(1)*(1.d0/3.d0 - dlnk/2.d0)))
  end function coop_cosmology_firstorder_sigma_Gaussian_R_quick
