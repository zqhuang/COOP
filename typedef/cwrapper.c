#ifdef HAS_FFTW
#include "fftw3.h"
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

void fft_1d_forward_(int *n,  double *in, fftw_complex *out){
  fftw_plan plan;
  plan = fftw_plan_dft_r2c_1d(*n, in, out, FFTW_ESTIMATE);  
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void fft_1d_backward_(int *n, fftw_complex *in, double *out){
  fftw_plan plan;
  plan = fftw_plan_dft_c2r_1d(*n, in, out, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void fft_2d_forward_(int *nx, int *ny,  double *in, fftw_complex *out){
  fftw_plan plan;
  plan = fftw_plan_dft_r2c_2d(*nx, *ny, in, out, FFTW_ESTIMATE);  //fortran switch
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void fft_2d_backward_(int *nx, int*ny, fftw_complex *in, double *out){
  fftw_plan plan;
  plan = fftw_plan_dft_c2r_2d(*nx, *ny, in, out, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void fft_3d_forward_(int *nx, int *ny, int *nz, double *in, fftw_complex *out){
  fftw_plan plan;
  plan = fftw_plan_dft_r2c_3d(*nx, *ny, *nz, in, out, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

void fft_3d_backward_(int *nx, int *ny, int *nz, fftw_complex *in, double *out){
  fftw_plan plan;
  plan = fftw_plan_dft_c2r_3d(*nx, *ny, *nz, in, out, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);}

#endif

void coop_array_copy_real_(double* x, double *y, int* n){
  int i;
  for(i=0; i<*n; i++)
    *(y+i) = *(x+i);
};

void coop_array_copy_single_(float* x, float *y, int* n){
  int i;
  for(i=0; i<*n; i++)
    *(y+i) = *(x+i);
};


void coop_array_copy_int_(int* x, int *y, int* n){
  int i;
  for(i=0; i<*n; i++)
    *(y+i) = *(x+i);
};


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
    if(fabs(*(x+i)) > mx){
      mx = fabs(*(x+i));
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
  if(i <0) i=0;
  if(buffer[i+1] > buffer[i]){
    r =  (*perc) + (1.*i)/(*n);
    *threshold = buffer[i-1]+(buffer[i]-buffer[i-1])*(1.-r);}
  else
    *threshold = buffer[i-1];
  free(buffer);
}

void array_get_mult_threshold_double_(double* x, int* n, double* perc, int* nthreshold, double* threshold){
  double *buffer;
  double r;
  int i, jt;
  buffer = (double*) malloc( (*n)*sizeof(double));
  for (i=0; i<*n; i++)
    buffer[i] = x[i];
  quicksort_double_(buffer, n);
  for(jt=0;jt< *nthreshold; jt++){
    i = (int) ((1. - (*(perc+jt)))*(*n)+0.5);
    if(i >= *n)i = *n-1;
    if(i <0) i=0;
    if(buffer[i+1] > buffer[i]){
      r =  (*(perc+jt)) + (1.*i)/(*n);
      *(threshold+jt) = buffer[i-1]+(buffer[i]-buffer[i-1])*(1.-r);}
    else
      *(threshold+jt)= buffer[i-1];}
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
  if(i <0) i=0;
  if(buffer[i+1] > buffer[i]){
    r =  (*perc) + (1.*i)/(*n);
    *threshold = buffer[i-1]+(buffer[i]-buffer[i-1])*(1.-r);}
  else
    *threshold = buffer[i-1];
  free(buffer);
}


void array_get_mult_threshold_float_(float* x, int* n, float* perc, int* nthreshold, float* threshold){
  float *buffer;
  float r;
  int i, jt;
  buffer = (float*) malloc( (*n)*sizeof(float));
  for (i=0; i<*n; i++)
    buffer[i] = x[i];
  quicksort_float_(buffer, n);
  for(jt=0;jt< *nthreshold; jt++){
    i = (int) ((1. - (*(perc+jt)))*(*n)+0.5);
    if(i >= *n)i = *n-1;
    if(i<0) i=0;
    if(buffer[i+1] > buffer[i]){
      r =  (*(perc+jt)) + (1.*i)/(*n);
      *(threshold+jt) = buffer[i-1]+(buffer[i]-buffer[i-1])*(1.-r);}
    else
      *(threshold+jt)= buffer[i-1];}
  free(buffer);
}


void cprint_(char *c){
  printf("%s", c);
}





void cwrapper_q_sortrev_int(int numbers[], int left, int right)
{
  int pivot, l_hold, r_hold;
#include "sortrev.h"
  if (l_hold < left)
    cwrapper_q_sortrev_int(numbers, l_hold, left-1);
  if (left < r_hold)
    cwrapper_q_sortrev_int(numbers, left+1, r_hold);
}

void quicksortrev_int_(int numbers[], int *array_size)
{
  cwrapper_q_sortrev_int(numbers, 0, (*array_size) - 1);
}


 
 
void cwrapper_q_sortrev_float(float numbers[], int left, int right)
{
  float pivot;
  int l_hold, r_hold;
#include "sortrev.h"
  if (l_hold < left)
    cwrapper_q_sortrev_float(numbers, l_hold, left-1);
  if (left < r_hold)
    cwrapper_q_sortrev_float(numbers, left+1, r_hold);
}


void quicksortrev_float_(float numbers[], int *array_size)
{
  cwrapper_q_sortrev_float(numbers, 0, (*array_size) - 1);
}

 
 
void cwrapper_q_sortrev_double(double numbers[], int left, int right)
{
  double pivot;
  int l_hold, r_hold;
#include "sortrev.h"
  if (l_hold < left)
    cwrapper_q_sortrev_double(numbers, l_hold, left-1);
  if (left < r_hold)
    cwrapper_q_sortrev_double(numbers, left+1, r_hold);
}


void quicksortrev_double_(double numbers[], int* array_size)
{
  cwrapper_q_sortrev_double(numbers, 0, (*array_size) - 1);
}



void cwrapper_q_sortrev_double_with_indices(double numbers[], int indices[], int left, int right)
{
  double pivot;
  int ipivot;
  int l_hold, r_hold;
#include "sortrev_acc.h"
  if (l_hold < left)
    cwrapper_q_sortrev_double_with_indices(numbers, indices, l_hold, left-1);
  if (left < r_hold)
    cwrapper_q_sortrev_double_with_indices(numbers, indices, left+1, r_hold);
}


void quicksortrev_double_with_indices_(double numbers[], int indices[], int* array_size)
{
  cwrapper_q_sortrev_double_with_indices(numbers, indices, 0, (*array_size) - 1);
}


void cwrapper_q_sortrev_float_with_indices(float numbers[], int indices[], int left, int right)
{
  float pivot;
  int ipivot;
  int l_hold, r_hold;
#include "sortrev_acc.h"
  if (l_hold < left)
    cwrapper_q_sortrev_float_with_indices(numbers, indices, l_hold, left-1);
  if (left < r_hold)
    cwrapper_q_sortrev_float_with_indices(numbers, indices, left+1, r_hold);
}


void quicksortrev_float_with_indices_(float numbers[], int indices[], int* array_size)
{
  cwrapper_q_sortrev_float_with_indices(numbers, indices, 0, (*array_size) - 1);
}


void cwrapper_q_sortrev_int_with_indices(int numbers[], int indices[], int left, int right)
{
  int pivot, ipivot;
  int l_hold, r_hold;
#include "sortrev_acc.h"
  if (l_hold < left)
    cwrapper_q_sortrev_int_with_indices(numbers, indices, l_hold, left-1);
  if (left < r_hold)
    cwrapper_q_sortrev_int_with_indices(numbers, indices, left+1, r_hold);
}


void quicksortrev_int_with_indices_(int numbers[], int indices[], int* array_size)
{
  cwrapper_q_sortrev_int_with_indices(numbers, indices, 0, (*array_size) - 1);
}




void cwrapper_q_sortrev_double_with_double(double numbers[], double indices[], int left, int right)
{
  double pivot;
  double ipivot;
  int l_hold, r_hold;
#include "sortrev_acc.h"
  if (l_hold < left)
    cwrapper_q_sortrev_double_with_double(numbers, indices, l_hold, left-1);
  if (left < r_hold)
    cwrapper_q_sortrev_double_with_double(numbers, indices, left+1, r_hold);
}


void quicksortrev_double_with_double_(double numbers[], double indices[], int* array_size)
{
  cwrapper_q_sortrev_double_with_double(numbers, indices, 0, (*array_size) - 1);
}


void cwrapper_q_sortrev_float_with_float(float numbers[], float indices[], int left, int right)
{
  float pivot;
  float ipivot;
  int l_hold, r_hold;
#include "sortrev_acc.h"
  if (l_hold < left)
    cwrapper_q_sortrev_float_with_float(numbers, indices, l_hold, left-1);
  if (left < r_hold)
    cwrapper_q_sortrev_float_with_float(numbers, indices, left+1, r_hold);
}


void quicksortrev_float_with_float_(float numbers[], float indices[], int* array_size)
{
  cwrapper_q_sortrev_float_with_float(numbers, indices, 0, (*array_size) - 1);
}
