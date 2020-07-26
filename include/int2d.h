    external func
    COOP_REAL integral
    COOP_REAL xmin, xmax, ymin, ymax, func, threshold1, threshold2
    COOP_INT, parameter::n_base = INT2D_N_BASE
    COOP_INT i, j
    COOP_REAL::int_grid(n_base, n_base), weight_grid(n_base, n_base), grids(0:n_base, 0:n_base), dx, dy, abs_int(n_base, n_base), center(n_base, n_base)
    weight_grid = 12.d0
    dx = (xmax-xmin)/n_base
    dy = (ymax-ymin)/n_base
    !$omp parallel do private(i, j)
    do i=0, n_base
       do j=0, n_base
          grids(i, j) =  func(xmin+dx*i, ymin+dy*j INT2D_ARGUMENTS)
       enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(i, j)
    do i=1, n_base
       do j=1, n_base
          center(i,j ) = func(xmin+(i-0.5d0)*dx, ymin+(j-0.5d0)*dy INT2D_ARGUMENTS)
          int_grid(i, j) =  (grids(i-1,j-1) + grids(i, j-1) + grids(i-1, j)+grids(i, j)) + 8.d0*center(i,j)
       enddo
    enddo
    !$omp end parallel do
    abs_int = abs(int_grid)
    call coop_array_get_threshold(abs_int, 0.01d0, threshold2)
    threshold1 = threshold2 * 0.2d0

    !$omp parallel do private(i, j)
    do i=1, n_base
       do j=1, n_base
          call refine_int(i, j)
       enddo
    enddo
    !$omp end parallel do

    integral = sum(int_grid/weight_grid)*dx*dy
    
  contains
    
    subroutine refine_int(i, j)
      COOP_INT i, j
      if(abs_int(i,j) .lt. threshold1)return
      if(abs_int(i,j) .lt. threshold2)then
         weight_grid(i,j) = weight_grid(i,j) + 36.d0 
         int_grid(i,j) = int_grid(i,j) - 4.d0*center(i,j) &
              + 2.d0*( &
              func(xmin+(i-1)*dx, ymin+(j-0.5d0)*dy INT2D_ARGUMENTS) &
              + func(xmin+i*dx, ymin+(j-0.5d0)*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-0.5d0)*dx, ymin+(j-1)*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-0.5d0)*dx, ymin+j*dy INT2D_ARGUMENTS) &
              ) &
              + 8.d0*( &
              func(xmin+(i-0.25d0)*dx, ymin + (j-0.25d0)*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-0.75d0)*dx, ymin + (j-0.25d0)*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-0.25d0)*dx, ymin + (j-0.75d0)*dy INT2D_ARGUMENTS)  &
              + func(xmin+(i-0.75d0)*dx, ymin + (j-0.75d0)*dy INT2D_ARGUMENTS)  &                              
              )
         return
      endif
      weight_grid(i,j) =  weight_grid(i,j) + 180.d0
      int_grid(i,j) = int_grid(i,j)  - 4.d0*center(i,j) &
           + 2.d0*( &
              func(xmin+(i-1)*dx, ymin+(j-0.5d0)*dy INT2D_ARGUMENTS) &
              + func(xmin+i*dx, ymin+(j-0.5d0)*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-0.5d0)*dx, ymin+(j-1)*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-0.5d0)*dx, ymin+j*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-1)*dx, ymin+(j-0.25d0)*dy INT2D_ARGUMENTS) &
              + func(xmin+i*dx, ymin+(j-0.25d0)*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-0.25d0)*dx, ymin+(j-1)*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-0.25d0)*dx, ymin+j*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-1)*dx, ymin+(j-0.75d0)*dy INT2D_ARGUMENTS) &
              + func(xmin+i*dx, ymin+(j-0.75d0)*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-0.75d0)*dx, ymin+(j-1)*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-0.75d0)*dx, ymin+j*dy INT2D_ARGUMENTS) &              
              ) + 4.d0*( &
              func(xmin+(i-0.25d0)*dx, ymin+(j-0.25d0)*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-0.5d0)*dx, ymin+(j-0.25d0)*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-0.75d0)*dx, ymin+(j-0.25d0)*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-0.25d0)*dx, ymin+(j-0.5d0)*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-0.75d0)*dx, ymin+(j-0.5d0)*dy INT2D_ARGUMENTS) &
              +func(xmin+(i-0.25d0)*dx, ymin+(j-0.75d0)*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-0.5d0)*dx, ymin+(j-0.75d0)*dy INT2D_ARGUMENTS) &
              + func(xmin+(i-0.75d0)*dx, ymin+(j-0.75d0)*dy INT2D_ARGUMENTS) &
              ) + 8.d0*( &
               func(xmin+(i-0.125d0)*dx, ymin+(j-0.125d0)*dy INT2D_ARGUMENTS) &
              +func(xmin+(i-0.375d0)*dx, ymin+(j-0.125d0)*dy INT2D_ARGUMENTS) &
              +func(xmin+(i-0.625d0)*dx, ymin+(j-0.125d0)*dy INT2D_ARGUMENTS) &
              +func(xmin+(i-0.875d0)*dx, ymin+(j-0.125d0)*dy INT2D_ARGUMENTS) &                             
              +func(xmin+(i-0.125d0)*dx, ymin+(j-0.375d0)*dy INT2D_ARGUMENTS) &
              +func(xmin+(i-0.375d0)*dx, ymin+(j-0.375d0)*dy INT2D_ARGUMENTS) &
              +func(xmin+(i-0.625d0)*dx, ymin+(j-0.375d0)*dy INT2D_ARGUMENTS) &
              +func(xmin+(i-0.875d0)*dx, ymin+(j-0.375d0)*dy INT2D_ARGUMENTS) &
              +func(xmin+(i-0.125d0)*dx, ymin+(j-0.625d0)*dy INT2D_ARGUMENTS) &
              +func(xmin+(i-0.375d0)*dx, ymin+(j-0.625d0)*dy INT2D_ARGUMENTS) &
              +func(xmin+(i-0.625d0)*dx, ymin+(j-0.625d0)*dy INT2D_ARGUMENTS) &
              +func(xmin+(i-0.875d0)*dx, ymin+(j-0.625d0)*dy INT2D_ARGUMENTS) &                             
              +func(xmin+(i-0.125d0)*dx, ymin+(j-0.875d0)*dy INT2D_ARGUMENTS) &
              +func(xmin+(i-0.375d0)*dx, ymin+(j-0.875d0)*dy INT2D_ARGUMENTS) &
              +func(xmin+(i-0.625d0)*dx, ymin+(j-0.875d0)*dy INT2D_ARGUMENTS) &
              +func(xmin+(i-0.875d0)*dx, ymin+(j-0.875d0)*dy INT2D_ARGUMENTS) &                                           
              )
      return
    end subroutine refine_int
