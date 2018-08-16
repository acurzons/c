#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "fitsio.h"

int main(){

  int arr = 8;
  int status, ii, jj, i;

  fitsfile *fptr;
  long fpixel, nelements;
  float *array[arr];
  
  char filename[50];
  char filefit[50] = "fitsfiles/rxj1713_Wp_ratio.fits";
  char filefit2[50];
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
  char paramsfile[200] = "rxj1713_Wp_ratio.txt";
  char paramsfile2[200];
  char buffer[250];
  char paramsline[100];
  float reg[29];

  params = fopen(paramsfile, "r");
  fgets(buffer, 200, params);

  for(i=0;i<29;i++){
    fgets(paramsline, sizeof(paramsline), params);
    sscanf(paramsline, "%f", &reg[i]);
   }

  fpixel = 1;
  nelements = (naxes[0]) * naxes[1];
  array[0] = (float *)malloc( naxes[0] * naxes[1] * sizeof(float));

  for(ii=0; ii<naxes[0]; ii++){
    array[ii] = array[0] + ii*naxes[1];
  }

  for(ii=0; ii<naxes[0]; ii++){
    for(jj=0; jj<naxes[1]; jj++){
      array[ii][jj] = 0.0;
    }
  }

  status = 0;

  remove(filefit);

  fits_create_file(&fptr, filefit, &status);

  fits_create_img(fptr, bitpix, naxis, naxes, &status);




  //Enter data here
  for(ii=0;ii<naxes[0];ii++){
    for(jj=0;jj<naxes[1];jj++){
      array[1][3] = reg[27];
      array[1][4] = reg[28];
      array[2][1] = reg[22];
      array[2][2] = reg[23];
      array[2][3] = reg[24];
      array[2][4] = reg[25];
      array[2][5] = reg[26];
      array[3][1] = reg[16];
      array[3][2] = reg[17];
      array[3][3] = reg[18];
      array[3][4] = reg[19];
      array[3][5] = reg[20];
      array[3][6] = reg[21];
      array[4][1] = reg[10];
      array[4][2] = reg[11];
      array[4][3] = reg[12];
      array[4][4] = reg[13];
      array[4][5] = reg[14];
      array[4][6] = reg[15];
      array[5][1] = reg[4];
      array[5][2] = reg[5];
      array[5][3] = reg[6];
      array[5][4] = reg[7];
      array[5][5] = reg[8];
      array[5][6] = reg[9];
      array[6][2] = reg[0];
      array[6][3] = reg[1];
      array[6][4] = reg[2];
      array[6][5] = reg[3];
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
