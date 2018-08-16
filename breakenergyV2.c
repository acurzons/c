#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "../newsedprod/src/constant.h"

int main(){

  int j;

  double B;
  double age = 5.05e10;
  double a = 637;
  double Ebreak;
  double Ecut;

    Ecut[j] = 230/pow(B[j],0.5);

    BG[j] = B[j] * 1e-6;

    Ebreak[j] = a/(pow(BG[j],2) * age);
    Ebreak[j] *= 0.62; //convert to TeV

    printf("Region %d has Eb = %.1lf, Ec = %.1lf and B = %.1lf \n", num, Ebreak[j], Ecut[j], B[j]);
 
  }
  
  fclose(fp);
  return(0);

}
