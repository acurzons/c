#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "fitsio.h"

int main(){

  int arr = 8;
  int status, ii, jj, kk, i;

  fitsfile *fptr;
  long fpixel, nelements;
  float *array[arr];
  
  char filefit[50];
  sprintf(filefit, "fitsfiles/rxj1713_sig_cor.fits");

  int bitpix = FLOAT_IMG;
  long naxis = 2;
  long naxes[2] = {arr, arr};

  char CTYPE1[100];
  sprintf(CTYPE1, "RA---CAR");
  char CTYPE2[100];
  sprintf(CTYPE2, "DEC--CAR");
  float CRPIX1 = 1;
  float CRPIX2 = 8;
  float CDELT1 = -0.18;
  float CDELT2 = 0.18;
  float CRVAL1 = 259.26;
  float CRVAL2 = -39.17;

  FILE *params;
  char paramsfile[200] = "rxj1713_regional_hadronic_parameters.txt";
  char buffer[250];
  char paramsline[100];
  float Wp[29];
  float buf;
  int reg[29];
  float outer, sig_neg, neg, no, pos, sig_pos;

  params = fopen(paramsfile, "r");
  fgets(buffer, 200, params);

  for(i=0;i<29;i++){
    fgets(paramsline, sizeof(paramsline), params);
    sscanf(paramsline, "%*d %*f %*f %f %*f %*f %*f %*f", &Wp[i]);
    //Wp[i] *= 1e46;
  }

  

  fpixel = 1;
  nelements = (naxes[0]) * naxes[1];
  array[0] = (float *)malloc( naxes[0] * naxes[1] * sizeof(float));

  for(ii=0; ii<naxes[0]; ii++){
    array[ii] = array[0] + ii*naxes[1];
  }


  outer = -4.0;
  for(ii=0; ii<naxes[0]; ii++){
    for(jj=0; jj<naxes[1]; jj++){
      array[ii][jj] = outer;
    }
  }

  status = 0;

  remove(filefit);

  fits_create_file(&fptr, filefit, &status);

  fits_create_img(fptr, bitpix, naxis, naxes, &status);

  sig_neg = -1.5;
  neg = -0.8;
  no = 1.1;
  pos = 1.5;
  sig_pos = 2.0;

  //Enter data here
  kk = 29;
  for(ii=0;ii<naxes[0];ii++){
    for(jj=0;jj<naxes[1];jj++){
      array[6][2] = sig_pos;
      array[6][3] = sig_neg;
      array[6][4] = no;
      array[6][5] = no;
      array[5][1] = pos;
      array[5][2] = no;
      array[5][3] = sig_pos;
      array[5][4] = no;
      array[5][5] = no;
      array[5][6] = sig_pos;
      array[4][1] = no;
      array[4][2] = no;
      array[4][3] = no;
      array[4][4] = no;
      array[4][5] = no;
      array[4][6] = no;
      array[3][1] = no;
      array[3][2] = no;
      array[3][3] = no;
      array[3][4] = no;
      array[3][5] = sig_pos;
      array[3][6] = no;
      array[2][1] = no;
      array[2][2] = no;
      array[2][3] = no;
      array[2][4] = no;
      array[2][5] = sig_pos;
      array[1][3] = no;
      array[1][4] = sig_pos;
      //  printf("array[ii][jj] = %.2f\n",array[ii][jj]);
    }
  }

  fits_write_img(fptr, TFLOAT, fpixel, nelements, array[0], &status);

  fits_update_key(fptr, TSTRING, "CTYPE1", CTYPE1, "Projection", &status);
  fits_update_key(fptr, TFLOAT, "CRPIX1", &CRPIX1, "reference pixel value", &status);
  fits_update_key(fptr, TFLOAT, "CDELT1", &CDELT1, "delta value per pixel", &status);
  fits_update_key(fptr, TFLOAT, "CRVAL1", &CRVAL1, "relative pixel value", &status);

  fits_update_key(fptr, TSTRING, "CTYPE2", CTYPE2, "Projection", &status);
  fits_update_key(fptr, TFLOAT, "CRPIX2", &CRPIX2, "reference pixel value", &status);
  fits_update_key(fptr, TFLOAT, "CDELT2", &CDELT2, "delta value per pixel", &status);
  fits_update_key(fptr, TFLOAT, "CRVAL2", &CRVAL2, "relative pixel value", &status);

  free(array[0]);
  
  fits_close_file(fptr, &status);
  

  return 0;
}
