#include <string.h>
#include <stdio.h>

#ifdef HAS_CFITSIO
#include "fitsio.h"
#endif

void coop_fits_read_header_to_string_(char* filename, char* cards, int* nkeys){
#ifdef HAS_CFITSIO
  fitsfile *fptr;         
  char card[FLEN_CARD]; 
  int status, ii; 
  status = 0; /* MUST initialize status */
  printf("%s", filename);
  fits_open_file(&fptr, filename, READONLY, &status);
  fits_get_hdrspace(fptr, nkeys, NULL, &status);
  for (ii = 1; ii <= *nkeys; ii++)  { 
    fits_read_record(fptr, ii, card, &status); /* read keyword */
    if(ii == 1)
      strcpy(cards, card);
    else
      strcat(cards, card);//    printf("%s\n", card);
    strcat(cards, "\n");//    printf("%s\n", card);
  }
  fits_close_file(fptr, &status);
    
  if (status)          /* print any error messages */
    fits_report_error(stderr, status);
#else
  printf("CFITSIO library is missing. Please change your configure file.");
#endif
}



void coop_fits_get_dimension_(char* filename, int* nx,int* ny, int* bytes){
#ifdef HAS_CFITSIO
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

#else
  printf("CFITSIO library is missing. Please change your configure file.");
#endif

}


void coop_fits_get_double_data_(char * filename, double* data, long * n){
#ifdef HAS_CFITSIO
  fitsfile *fptr;         
  int status;
  long fpixel[2];
  LONGLONG nelements;
   fpixel[0] = 1;
   fpixel[1] = 1;
   status = 0; /* MUST initialize status */
   nelements = *n;
   fits_open_file(&fptr, filename, READONLY, &status);
   fits_read_pix(fptr, TDOUBLE, fpixel, nelements, NULL, data, NULL, &status);  
   fits_close_file(fptr, &status);
   if (status)          /* print any error messages */
     fits_report_error(stderr, status);
#else
  printf("CFITSIO library is missing. Please change your configure file.");
#endif
}


void coop_fits_get_float_data_(char * filename, float* data, long * n){
#ifdef HAS_CFITSIO
  fitsfile *fptr;         
  int status;
  long fpixel[2];
  LONGLONG nelements;
   fpixel[0] = 1;
   fpixel[1] = 1;
   status = 0; /* MUST initialize status */
   nelements = *n;
   fits_open_file(&fptr, filename, READONLY, &status);
   fits_read_pix(fptr, TDOUBLE, fpixel, nelements, NULL, data, NULL, &status);  
   fits_close_file(fptr, &status);
   if (status)          /* print any error messages */
     fits_report_error(stderr, status);
#else
  printf("CFITSIO library is missing. Please change your configure file.");
#endif
}

