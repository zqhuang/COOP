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
#include <time.h>
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


 
void cwrapper_q_sort_int(int numbers[], int left, int right)
{
  int pivot, l_hold, r_hold;
#include "sort.h"
  if (l_hold < left)
    cwrapper_q_sort_int(numbers, l_hold, left-1);
  if (left < r_hold)
    cwrapper_q_sort_int(numbers, left+1, r_hold);
}

void quicksort_int_(int numbers[], int *array_size)
{
  cwrapper_q_sort_int(numbers, 0, (*array_size) - 1);
}


 
 
void cwrapper_q_sort_float(float numbers[], int left, int right)
{
  float pivot;
  int l_hold, r_hold;
#include "sort.h"
  if (l_hold < left)
    cwrapper_q_sort_float(numbers, l_hold, left-1);
  if (left < r_hold)
    cwrapper_q_sort_float(numbers, left+1, r_hold);
}


void quicksort_float_(float numbers[], int *array_size)
{
  cwrapper_q_sort_float(numbers, 0, (*array_size) - 1);
}

 
 
void cwrapper_q_sort_double(double numbers[], int left, int right)
{
  double pivot;
  int l_hold, r_hold;
#include "sort.h"
  if (l_hold < left)
    cwrapper_q_sort_double(numbers, l_hold, left-1);
  if (left < r_hold)
    cwrapper_q_sort_double(numbers, left+1, r_hold);
}


void quicksort_double_(double numbers[], int* array_size)
{
  cwrapper_q_sort_double(numbers, 0, (*array_size) - 1);
}



void cwrapper_q_sort_double_with_indices(double numbers[], int indices[], int left, int right)
{
  double pivot;
  int ipivot;
  int l_hold, r_hold;
#include "sort_acc.h"
  if (l_hold < left)
    cwrapper_q_sort_double_with_indices(numbers, indices, l_hold, left-1);
  if (left < r_hold)
    cwrapper_q_sort_double_with_indices(numbers, indices, left+1, r_hold);
}


void quicksort_double_with_indices_(double numbers[], int indices[], int* array_size)
{
  cwrapper_q_sort_double_with_indices(numbers, indices, 0, (*array_size) - 1);
}


void cwrapper_q_sort_float_with_indices(float numbers[], int indices[], int left, int right)
{
  float pivot;
  int ipivot;
  int l_hold, r_hold;
#include "sort_acc.h"
  if (l_hold < left)
    cwrapper_q_sort_float_with_indices(numbers, indices, l_hold, left-1);
  if (left < r_hold)
    cwrapper_q_sort_float_with_indices(numbers, indices, left+1, r_hold);
}


void quicksort_float_with_indices_(float numbers[], int indices[], int* array_size)
{
  cwrapper_q_sort_float_with_indices(numbers, indices, 0, (*array_size) - 1);
}


void cwrapper_q_sort_int_with_indices(int numbers[], int indices[], int left, int right)
{
  int pivot, ipivot;
  int l_hold, r_hold;
#include "sort_acc.h"
  if (l_hold < left)
    cwrapper_q_sort_int_with_indices(numbers, indices, l_hold, left-1);
  if (left < r_hold)
    cwrapper_q_sort_int_with_indices(numbers, indices, left+1, r_hold);
}


void quicksort_int_with_indices_(int numbers[], int indices[], int* array_size)
{
  cwrapper_q_sort_int_with_indices(numbers, indices, 0, (*array_size) - 1);
}




void cwrapper_q_sort_double_with_double(double numbers[], double indices[], int left, int right)
{
  double pivot;
  double ipivot;
  int l_hold, r_hold;
#include "sort_acc.h"
  if (l_hold < left)
    cwrapper_q_sort_double_with_double(numbers, indices, l_hold, left-1);
  if (left < r_hold)
    cwrapper_q_sort_double_with_double(numbers, indices, left+1, r_hold);
}


void quicksort_double_with_double_(double numbers[], double indices[], int* array_size)
{
  cwrapper_q_sort_double_with_double(numbers, indices, 0, (*array_size) - 1);
}


void cwrapper_q_sort_float_with_float(float numbers[], float indices[], int left, int right)
{
  float pivot;
  float ipivot;
  int l_hold, r_hold;
#include "sort_acc.h"
  if (l_hold < left)
    cwrapper_q_sort_float_with_float(numbers, indices, l_hold, left-1);
  if (left < r_hold)
    cwrapper_q_sort_float_with_float(numbers, indices, left+1, r_hold);
}


void quicksort_float_with_float_(float numbers[], float indices[], int* array_size)
{
  cwrapper_q_sort_float_with_float(numbers, indices, 0, (*array_size) - 1);
}

void get_binned_data_double_(double x[], int* n, int* nbins, int* bins_used, double center[], double density[]){
  double* buffer;
  double lower, upper, dx;
  int i, m, j, i1;
  if(*n < 2){
    center[0]= x[0];
    *bins_used = 1;
    density[0] = 1.;
    return;}
  buffer = (double*) malloc((*n)*sizeof(double));
  for (i=0; i<*n; i++)
    buffer[i] = x[i];
  quicksort_double_(buffer, n);
  m = (*n - 1)/(*nbins) + 1;
  lower = buffer[0]-1.e-30;
  upper = lower;
  dx = (buffer[*n-1] - buffer[0])/(*nbins)/2.;
  i = 0;
  j = -1;
  i1 = i;
  while(i < *n){
    j++;
    center[j] = 0;
    while(upper < lower+dx  || i1-i <m){
      center[j]= center[j]+ buffer[i1];
      if(i1 == (*n) - 1){
	upper = buffer[i1]+1.e-30;
	i1++ ;
	break;}
      else{
	upper = (buffer[i1]+buffer[i1+1])/2.;
	i1++;}
    }
    center[j] = center[j]/(i1 - i);
    density[j] = (i1-i)/(upper - lower);
    lower = upper;
    i = i1;
  }
  free(buffer);
  *bins_used = j+1;
  if(density[0] > density[1])
    density[0] = density[1];
  if(density[j] > density[j-1])
    density[j] = density[j-1];
}


void get_binned_data_float_(float x[], int* n, int* nbins, int* bins_used, float center[], float density[]){
  float* buffer;
  float lower, upper, dx;
  int i, m, j, i1;
  if(*n < 2){
    center[0]= x[0];
    *bins_used = 1;
    density[0] = 1.;
    return;}
  buffer = (float*) malloc((*n)*sizeof(float));
  for (i=0; i<*n; i++)
    buffer[i] = x[i];
  quicksort_float_(buffer, n);
  m = (*n - 1)/(*nbins) + 1;
  lower = buffer[0]-1.e-30;
  upper = lower;
  dx = (buffer[*n-1] - buffer[0])/(*nbins)/2.;
  i = 0;
  j = -1;
  i1 = i;
  while(i < *n){
    j++;
    center[j] = 0;
    while(upper < lower+dx  || i1-i <m){
      center[j]= center[j]+ buffer[i1];
      if(i1 == (*n) - 1){
	upper = buffer[i1]+1.e-30;
	i1++ ;
	break;}
      else{
	upper = (buffer[i1]+buffer[i1+1])/2.;
	i1++;}
    }
    center[j] = center[j]/(i1 - i);
    density[j] = (i1-i)/(upper - lower);
    lower = upper;
    i = i1;
  }
  free(buffer);
  *bins_used = j+1;
  if(density[0] > density[1])
    density[0] = density[1];
  if(density[j] > density[j-1])
    density[j] = density[j-1];
}


//find threshold such that P(x>threashold) is perc
void array_get_threshold_double_(double* x, int* n, double* perc, double* threshold){
  double *buffer;
  double r;
  int i;
  buffer = (double*) malloc( (*n)*sizeof(double));
  for (i=0; i<*n; i++)
    buffer[i] = x[i];
  quicksort_double_(buffer, n);
  i = (int) ((1. - (*perc))*(*n)+0.5);
  if(i >= *n)i = *n-1;
  if(buffer[i+1] > buffer[i]){
    r =  (*perc) + (1.*i)/(*n);
    *threshold = buffer[i-1]+(buffer[i]-buffer[i-1])*(1.-r);}
  else
    *threshold = buffer[i-1];
  free(buffer);
}

//find threshold such that P(x>threashold) is perc
void array_get_threshold_float_(float* x, int* n, float* perc, float* threshold){
  float *buffer;
  float r;
  int i;
  buffer = (float*) malloc( (*n)*sizeof(float));
  for (i=0; i<*n; i++)
    buffer[i] = x[i];
  quicksort_float_(buffer, n);
  i = (int) ((1. - (*perc))*(*n)+0.5);
  if(i >= *n)i = *n-1;
  if(buffer[i+1] > buffer[i]){
    r =  (*perc) + (1.*i)/(*n);
    *threshold = buffer[i-1]+(buffer[i]-buffer[i-1])*(1.-r);}
  else
    *threshold = buffer[i-1];
  free(buffer);
}




 
