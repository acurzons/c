#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "../newsedprod/src/constant.h"

int main(){

  FILE *fp;
  char* filename = "rxj1713_leptonic_regional_parameters.txt";
  char buffer[250];

  int j;

  double mass = 9.11e-31;
  double e = 1.6e-19;
  double c = 3.0e8;

  int num;
  double N_H;
  double n_H[29];
  double W[29];
  double W1TeV;
  double B[29];
  double a1[29];
  double a2[29];
  double Ec[29];
  double Eb;
  double BG[29];

  double age = 5.05e10;
  double a = 637;
  double Ebreak[29];
  double Ecut[29];

  fp = fopen(filename, "r");
  fgets(buffer, 150, fp);

  for(j=0;j<29;j++){
    fscanf(fp, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",&num,&N_H,&n_H[j],&W[j],&W1TeV,&B[j],&a1[j],&a2[j],&Ec[j],&Eb);

    Ecut[j] = 230/pow(B[j],0.5);

    BG[j] = B[j] * 1e-6;

    Ebreak[j] = a/(pow(BG[j],2) * age);
    Ebreak[j] *= 0.62; //convert to TeV

    printf("Region %d has Eb = %.1lf, Ec = %.1lf and B = %.1lf \n", num, Ebreak[j], Ecut[j], B[j]);
 
  }
  
  fclose(fp);
  return(0);

}
