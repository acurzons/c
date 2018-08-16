#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "../newsedprod/src/constant.h"


int main(){

  FILE *params;
  char paramsfilename[250];
  char paramsline[100];

  int j;
 
  double E;
  double flux;
  double cumsum = 0.0;
 
  for(j=0;j<29;j++){

    sprintf(paramsfilename, "../newsedprod/OUTPUT/rxj1713_had%d_ppgamma.txt", j+1);
    params = fopen(paramsfilename, "r");

    while(fgets(paramsline, sizeof(paramsline), params)){
      sscanf(paramsline, "%lf %*lf %*lf %*lf %*lf %lf", &E, &flux);
      if(E == 1.259e-3) cumsum += flux;
    }

    
    
    fclose(params);
  }  
  
  printf("Cumsum = %.3e\n", cumsum);
  
return(0);

}
