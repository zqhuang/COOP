#include <string.h>
#include <stdio.h>
#include "fitsio.h"

void coop_fits_print_header_(char* filename){
  fitsfile *fptr;         
  char card[FLEN_CARD]; 
  int status,  nkeys, ii;  
  status = 0; /* MUST initialize status */
  fits_open_file(&fptr, filename, READONLY, &status);
  fits_get_hdrspace(fptr, &nkeys, NULL, &status);
  
  for (ii = 1; ii <= nkeys; ii++)  { 
    fits_read_record(fptr, ii, card, &status); /* read keyword */
    printf("%s\n", card);
  }
  printf("END\n\n");  /* terminate listing with END */
  fits_close_file(fptr, &status);
  
  if (status)          /* print any error messages */
    fits_report_error(stderr, status);
}


void coop_fits_get_dimension_(char* filename, int* nx,int* ny, int* bytes){
  fitsfile *fptr;         
  int status, idim;
  long naxes[2];
  status = 0; /* MUST initialize status */
  fits_open_file(&fptr, filename, READONLY, &status);
  fits_get_img_dim(fptr, & idim, & status);
  if(idim != 2) {
    *nx = 0;
    *ny = 0;
    *bytes = 0;
  }
  else{
    fits_get_img_size(fptr, idim, naxes, &status);
    *nx = naxes[0];
    *ny = naxes[1];
    fits_get_img_type(fptr, bytes, &status);
    *bytes = - (*bytes) / 8;
  }
  fits_close_file(fptr, &status);
  if (status)          /* print any error messages */
    fits_report_error(stderr, status);

}

void coop_fits_get_float_data_(char * filename, float * data, int * nx, int* ny){
  fitsfile *fptr;         
  int status, idim;
  long fpixel[2];
  long nelements;
  float nval;
  int * anynul;
  fpixel[0] = 0;
  fpixel[1] = 0;
  status = 0; /* MUST initialize status */
  nval = 0.;
  nelements = (*nx)*(*ny);
  fits_open_file(&fptr, filename, READONLY, &status);
  fits_read_pix(fptr, TFLOAT, fpixel,nelements, &nval, data, anynul, &status);  
  fits_close_file(fptr, &status);
  if (status)          /* print any error messages */
    fits_report_error(stderr, status);
}

void coop_fits_get_double_data_(char * filename, double* data, int * nx, int* ny){
  fitsfile *fptr;         
  int status;
  long fpixel[2];
  long nelements;
  double nval;
  int* anynul;
   fpixel[0] = 1;
   fpixel[1] = 1;
   status = 0; /* MUST initialize status */
    nval = 0.;
   nelements = (*nx)*(*ny);
   fits_open_file(&fptr, filename, READONLY, &status);
     fits_read_pix(fptr, TDOUBLE, fpixel, nelements, &nval, data, anynul, &status);  
   fits_close_file(fptr, &status);
   if (status)          /* print any error messages */
     fits_report_error(stderr, status);
}




/*int main(int argc, char *argv[])
{
  int nx, ny;
  fits_get_dimension_(argv[1], &nx, &ny);
  return nx;
 }
*/
