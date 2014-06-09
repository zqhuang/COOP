#ifdef HAS_FFTW
#include "fftw3.h"
#ifdef MPI
#include "fftw3-mpi.h"
#endif
#endif
#ifdef MPI
#include "mpi.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef HAS_GSL
#include "gsl_interp.h"
#include "gsl_spline.h"
#include "gsl_math.h"
#include "gsl_sf.h"
#endif
#ifdef HAS_FFTW
//********************* fft ***************************************//
void fft_dcti_(int* n, double* in, double* out){
  fftw_plan plan;
  plan = fftw_plan_r2r_1d(*n, in, out, FFTW_REDFT00, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void fft_dsti_(int* n, double* in, double* out){
  fftw_plan plan;
  plan = fftw_plan_r2r_1d(*n, in, out, FFTW_RODFT00, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void fft_dctii_(int* n, double* in, double* out){
  fftw_plan plan;
  plan = fftw_plan_r2r_1d(*n, in, out, FFTW_REDFT10, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void fft_dstii_(int* n, double* in, double* out){
  fftw_plan plan;
  plan = fftw_plan_r2r_1d(*n, in, out, FFTW_RODFT10, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void fft_idctii_(int* n, double* in, double* out){
  fftw_plan plan;
  plan = fftw_plan_r2r_1d(*n, in, out, FFTW_REDFT01, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
}

void fft_idstii_(int* n, double* in, double* out){
  fftw_plan plan;
  plan = fftw_plan_r2r_1d(*n, in, out, FFTW_RODFT01, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}


void fft_1d_forward_(int *n0, double *in, fftw_complex *out){
  fftw_plan plan;
  plan = fftw_plan_dft_r2c_1d(*n0, in, out, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void fft_1d_backward_(int *n0, fftw_complex *in, double *out){
  fftw_plan plan;
  plan = fftw_plan_dft_c2r_1d(*n0, in, out, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}


void fft_2d_forward_(int *n0, double *in, fftw_complex *out){
  fftw_plan plan;
  plan = fftw_plan_dft_r2c_2d(*n0, *n0, in, out, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void fft_2d_backward_(int *n0, fftw_complex *in, double *out){
  fftw_plan plan;
  plan = fftw_plan_dft_c2r_2d(*n0, *n0, in, out, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void fft_3d_forward_(int *n0, double *in, fftw_complex *out){
  fftw_plan plan;
  plan = fftw_plan_dft_r2c_3d(*n0, *n0, *n0, in, out, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void fft_3d_backward_(int *n0, fftw_complex *in, double *out){
  fftw_plan plan;
  plan = fftw_plan_dft_c2r_3d(*n0, *n0, *n0, in, out, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void inplacefft_3d_forward_(int *n0, double *in){
  fftw_plan plan;
  plan = fftw_plan_dft_r2c_3d(*n0, *n0, *n0, in, (fftw_complex *) in, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void inplacefft_3d_backward_(int *n0, double *in){
  fftw_plan plan;
  plan = fftw_plan_dft_c2r_3d(*n0, *n0, *n0, (fftw_complex *) in,  in, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void inplacefft_3d_forward_saved_(int *n0, double *in){
  static fftw_plan plan;
  if(*n0 >0)
    plan = fftw_plan_dft_r2c_3d(*n0, *n0, *n0, in, (fftw_complex *) in, FFTW_MEASURE);
  else if(*n0 == 0)
    fftw_execute(plan);
  else
    fftw_destroy_plan(plan);}

void inplacefft_3d_backward_saved_(int *n0, double *in){
  static fftw_plan plan ;
  if(*n0 > 0)
    plan = fftw_plan_dft_c2r_3d(*n0, *n0, *n0, (fftw_complex *) in,  in, FFTW_MEASURE);
  else if(*n0 ==0)
    fftw_execute(plan);
  else
    fftw_destroy_plan(plan);}


#ifdef MPI

void mpifft_init_(void){
  fftw_mpi_init();}

//the size of the input array should be n0 x n0 x [2(n0/2 + 1)] globally
//locally the input array would be n0_local x n0 x [2(n0/2 +1)]
void inplacempifft_3d_forward_(int *n, double *in){
  ptrdiff_t n0;
  fftw_plan plan;
  n0 = *n;
  plan = fftw_mpi_plan_dft_r2c_3d(n0, n0, n0, in, (fftw_complex *) in, MPI_COMM_WORLD, FFTW_MPI_TRANSPOSED_OUT | FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

//the size of the input array should be n0 x n0 x [2(n0/2 + 1)] globally
//locally the input array would be n0_local x n0 x [2(n0/2 +1)] 
void inplacempifft_3d_backward_(int *n, double *in){
  fftw_plan plan;
  ptrdiff_t n0;
  n0 = *n;
  plan = fftw_mpi_plan_dft_c2r_3d(n0, n0, n0, (fftw_complex *) in, in, MPI_COMM_WORLD, FFTW_MPI_TRANSPOSED_IN | FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void inplacempifft_3d_forward_saved_(int *n, double *in){
  ptrdiff_t n0;
  static fftw_plan plan;
  if(*n >0){
    n0 = *n;
    plan = fftw_mpi_plan_dft_r2c_3d(n0, n0, n0, in, (fftw_complex *) in, MPI_COMM_WORLD, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);}
  else if(*n == 0)
    fftw_execute(plan);
  else
    fftw_destroy_plan(plan);}

//the size of the input array should be n0 x n0 x [2(n0/2 + 1)] globally
//locally the input array would be n0_local x n0 x [2(n0/2 +1)] 
void inplacempifft_3d_backward_saved_(int *n, double *in){
  static fftw_plan plan;
  ptrdiff_t n0;
  if(*n>0){
    n0 = *n;
    plan = fftw_mpi_plan_dft_c2r_3d(n0, n0, n0, (fftw_complex *) in, in, MPI_COMM_WORLD, FFTW_MPI_TRANSPOSED_IN | FFTW_MEASURE);}
  else if( *n ==0)
    fftw_execute(plan);
  else
    fftw_destroy_plan(plan);}



void mpifft_transpose_get_local_n_(int *n_total, int *n_local_x, int *n_local_x_start, int *n_local_y, int *n_local_y_start, int *nz_local){
  ptrdiff_t n0, local_n0, local_0_start, local_n1, local_1_start, alloc_local;
  n0 = *n_total;
  alloc_local = fftw_mpi_local_size_3d_transposed(n0, n0, n0/2+1, MPI_COMM_WORLD, &local_n0, &local_0_start, &local_n1, &local_1_start);
  *nz_local = alloc_local/((n0/2+1) * n0);
  if(alloc_local > (*nz_local) * ((n0/2+1) * n0) ) (*nz_local)++;
  *n_local_x = local_n0;
  *n_local_x_start = local_0_start+1; //Fortran index
  *n_local_y = local_n1;
  *n_local_y_start = local_1_start+1; //Fortran index}  
}
#endif


#endif



#ifdef HAS_GSL
void spherical_bessel_jl_(int* l, double* x, double* result){
  *result = gsl_sf_bessel_jl(*l, *x);
};

void spherical_bessel_jl_array_(int *lmax, double *x, double *result_array){
  gsl_sf_bessel_jl_steed_array(*lmax, *x, result_array);
};

void spherical_harmonics_ylm_(int *l, int* m, double* x, double* Ylm){
  *Ylm = gsl_sf_legendre_sphPlm(*l, *m, *x);
  printf("%d %d %e %e", *l,*m, *x, *Ylm);
};
#endif


void count_array_threshold_(double* x, int* n, double* threshold, int* nlarge){
  int i;
  *nlarge = 0;
  for(i=0; i<*n; i++)
    if(*(x+i) > *threshold) (*nlarge)++;
};


void array_maxabs_location_(double* x, int* n, int* loc){
  int i;
  double mx;
  mx = 0;
  for(i=0;i< *n;i++){
    if(abs(*(x+i)) > mx){
      mx = abs(*(x+i));
      *loc = i;
    }
  }
}

void array_get_maxval_(double* x, int* n, double *xmax){
  int i;
  *xmax = x[0];
  for(i=1;i< *n;i++){
    if( x[i] > *xmax)
      *xmax = x[i];
  }
}


void array_get_minval_(double* x, int* n, double *xmin){
  int i;
  *xmin = x[0];
  for(i=1;i< *n;i++){
    if(x[i] < *xmin)
      *xmin = x[i];
  }
}

void float_array_get_maxval_(float* x, int* n, float *xmax){
  int i;
  *xmax = x[0];
  for(i=1;i< *n;i++){
    if( x[i] > *xmax)
      *xmax = x[i];
  }
}


void float_array_get_minval_(float* x, int* n, float *xmin){
  int i;
  *xmin = x[0];
  for(i=1;i< *n;i++){
    if(x[i] < *xmin)
      *xmin = x[i];
  }
}


//find threshold such that P(x>threashold) is perc
void array_get_threshold_(double* x, int* n, double* perc, double* threshold){
  const int m = 562144;
  int i, ncut, nacc;
  double xmin, xmax, step;
  int a[m];
  array_get_maxval_(x, n, &xmax);
  array_get_minval_(x, n, &xmin);
  if(*perc > 1. || xmax==xmin){
    *threshold = xmin;
    return;
  }
  if(*perc < 0.){
    *threshold = xmax;
    return;
  }
  step = (xmax-xmin)/(m-1); 
  for(i=0; i<m; i++)
    a[i] = 0;
  for(i=0; i< *n; i++)
    (a[(int)((x[i] - xmin)/step+0.5)])++; 
  ncut = (int) ((*n)*(*perc)+0.5);
  i = m-1;
  nacc = a[i]; 
  while( nacc < ncut)
    nacc += a[--i];
  if(i == 0)
    *threshold = xmin + 0.25*step;
  else if(i==m-1)
    *threshold = xmax - 0.25*step;
  else
    *threshold = i*step+xmin;
}


void float_array_get_threshold_(float* x, int* n, float* perc, float* threshold){
  const int m = 262144;
  int i, ncut, nacc;
  float xmin, xmax, step;
  int a[m];
  float_array_get_maxval_(x, n, &xmax);
  float_array_get_minval_(x, n, &xmin);
  if(*perc > 1. || xmax==xmin){
    *threshold = xmin;
    return;
  }
  if(*perc < 0.){
    *threshold = xmax;
    return;
  }
  step = (xmax-xmin)/(m-1); 
  for(i=0; i<m; i++)
    a[i] = 0;
  for(i=0; i< *n; i++)
    (a[(int)((x[i] - xmin)/step+0.5)])++; 
  ncut = (int) ((*n)*(*perc)+0.5);
  i = m-1;
  nacc = a[i]; 
  while( nacc < ncut)
    nacc += a[--i];
  if(i == 0)
    *threshold = xmin + 0.25*step;
  else if(i==m-1)
    *threshold = xmax - 0.25*step;
  else
    *threshold = i*step+xmin;
}


 
 
void cwrapper_q_sort_int(int numbers[], int left, int right)
{
  int pivot, l_hold, r_hold;
  l_hold = left;
  r_hold = right;
  pivot = numbers[left];
  while (left < right)
  {
    while ((numbers[right] >= pivot) && (left < right))
      right--;
    if (left != right)
    {
      numbers[left] = numbers[right];
      left++;
    }
    while ((numbers[left] <= pivot) && (left < right))
      left++;
    if (left != right)
    {
      numbers[right] = numbers[left];
      right--;
    }
  }
  numbers[left] = pivot;
  pivot = left;
  left = l_hold;
  right = r_hold;
  if (left < pivot)
    cwrapper_q_sort_int(numbers, left, pivot-1);
  if (right > pivot)
    cwrapper_q_sort_int(numbers, pivot+1, right);
}

void quicksort_int_(int numbers[], int *array_size)
{
  cwrapper_q_sort_int(numbers, 0, (*array_size) - 1);
}


 
 
void cwrapper_q_sort_float(float numbers[], int left, int right)
{
  float pivot;
  int l_hold, r_hold;
 
  l_hold = left;
  r_hold = right;
  pivot = numbers[left];
  while (left < right)
  {
    while ((numbers[right] >= pivot) && (left < right))
      right--;
    if (left != right)
    {
      numbers[left] = numbers[right];
      left++;
    }
    while ((numbers[left] <= pivot) && (left < right))
      left++;
    if (left != right)
    {
      numbers[right] = numbers[left];
      right--;
    }
  }
  numbers[left] = pivot;
  right = r_hold;
  r_hold = left;
  left = l_hold;
  if (left < r_hold)
    cwrapper_q_sort_float(numbers, left, r_hold-1);
  if (right > r_hold)
    cwrapper_q_sort_float(numbers, r_hold+1, right);
}


void quicksort_float_(float numbers[], int *array_size)
{
  cwrapper_q_sort_float(numbers, 0, (*array_size) - 1);
}

 
 
void cwrapper_q_sort_double(double numbers[], int left, int right)
{
  double pivot;
  int l_hold, r_hold;
 
  l_hold = left;
  r_hold = right;
  pivot = numbers[left];
  while (left < right)
  {
    while ((numbers[right] >= pivot) && (left < right))
      right--;
    if (left != right)
    {
      numbers[left] = numbers[right];
      left++;
    }
    while ((numbers[left] <= pivot) && (left < right))
      left++;
    if (left != right)
    {
      numbers[right] = numbers[left];
      right--;
    }
  }
  numbers[left] = pivot;
  right = r_hold;
  r_hold = left;
  left = l_hold;
  if (left < r_hold)
    cwrapper_q_sort_double(numbers, left, r_hold-1);
  if (right > r_hold)
    cwrapper_q_sort_double(numbers, r_hold+1, right);
}


void quicksort_double_(double numbers[], int* array_size)
{
  cwrapper_q_sort_double(numbers, 0, (*array_size) - 1);
}
