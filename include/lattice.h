#define DO_LOOP  do k=0, this%n-1; do j=0,this%n-1; do i=0, this%n-1
#define END_LOOP  enddo;enddo;enddo
#define LATTICE_INDS i, j, k
#define F_ALL  this%f(:, 0:this%n-1, 0:this%n-1, 0:this%n-1)
#define V_INDEX  (2.d0*(3.d0+this%mu)/(3.d0-this%mu))
#define G_INDEX  (2.d0*(1.d0+this%mu)/(3.d0-this%mu))
#define K_INDEX  (-2.d0)
#define A_INDEX  (2.d0/(3.d0-this%mu))
#define GR_COEF (-(3.d0-this%mu)**2/(24.d0*coop_lattice_mpsq))
#define LAP_F(i, j, k) (this%f(:, i+1, j, k) + this%f(:, i-1, j, k) + this%f(:, i, j+1, k) + this%f(:, i, j-1, k) + this%f(:, i, j, k-1) + this%f(:, i, j, k+1) - 6.d0*this%f(:,i, j, k))

#define IND11 1
#define IND22 2
#define IND33 3
#define IND23 4
#define IND32 4
#define IND31 5
#define IND13 5
#define IND12 6
#define IND21 6
