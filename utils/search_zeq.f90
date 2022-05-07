program searchzeq
  use parabolaU_mod
  double precision :: z_eq, z2
  integer::i
  do i=1, 10
     z2 = 5.1d0+i*0.01d0
     z_eq = z_equilibrium(focal_length=2.1d0, z_gap12 = 0.1d0, z_gap23=z2, z_max=29.195d0, size_gap12=0.1d0, size_gap23=0.1d0, potential_1=0.d0, potential_2=-1.d0, potential_3=0.d0)
     write(*,*) z2, z_eq
  enddo
  
end program searchzeq

