#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "../newsedprod/src/constant.h"

int main(){

  double B = 0.0;
  double BG;
  double age = 5.05e10;
  double a = 637;
  double Ebreak;
  double Ecut;

  printf("Please enter B in units of microGauss\n");
  scanf("%lf",&B);

  Ecut = 230/pow(B,0.5);

  BG = B * 1e-6;

  Ebreak = a/(pow(BG,2) * age);
  Ebreak *= 0.62; //convert to TeV

  printf("Eb = %.1lf TeV\nEc = %.1lf TeV\nB = %.1lf microGauss \n", Ebreak, Ecut, B);
  
  return(0);

}
