#define LATTICE_INDS  ii, jj, kk, kmin, kmax, jmin, jmax, imin, imax, i, j, k, fld
#define  LOOP_LEVEL1  do kk=0, this%nc-1; kmin=kk*this%cs; kmax = kmin + this%cs -1; do jj = 0, this%nc-1; jmin = jj*this%cs; jmax = jmin + this%cs-1; do ii = 0, this%nc-1; imin = ii*this%cs; imax = imin + this%cs-1
#define LOOP_LEVEL2  do k=kmin, kmax; do j=jmin, jmax; do i = imin, imax
#define END_LOOP enddo; enddo; enddo

#define IND11 1
#define IND22 2
#define IND33 3
#define IND23 4
#define IND32 4
#define IND31 5
#define IND13 5
#define IND12 6
#define IND21 6
