#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"

int main(){

  FILE *brems;
  FILE *IC;
  
  double brems_flux[561];
  double brems_energy[561];
  double brems_extra[10000];
  double IC_flux[561];
  double IC_energy[561];
  double IC_extra[10000];

  int i = 0;
  int size = 0;

  IC = fopen("../newsedprod/OUTPUT/rxj1713_lep1_bpl_ICprimary.txt","r");
  brems = fopen("../newsedprod/OUTPUT/rxj1713_lep1_bpl_bremssprimary.txt","r");
  
  while(!feof(IC)){
    fscanf(IC, "%lf", &IC_extra[size]);
    size++;
  }
  printf("size = %d\n",size);

  /*
  fseek(brems, 0, SEEK_END);
  size = ftell(brems);
  fseek(brems, 0, SEEK_SET);
  size /= sizeof(double);
  printf("size = %d\n",size);
  */
 

  for(i=0; i<size2; i++){
    fscanf(brems, "%lf %*lf %*lf %*lf %*lf %lf", &brems_energy[i], &brems_flux[i]);
    printf("i = %d, flux = %.3e\n",i, brems_flux[i/6]);
    
    /*    if(i%6 == 0){
      fscanf(brems, "%lf", &brems_flux[i/6]);
      printf("i = %d, flux = %.3e\n",i, brems_flux[i/6]);
    }
    else if((i+1)%6 == 0){
      fscanf(brems, "%lf", &brems_energy[(i - 5)/6]);
  }
    else{
      fscanf(brems, "%lf", &brems_extra[i]);
      }*/
  }
  printf("i = %d\n",i);

  /*
  i = 0;
  for(i=0; i<size2; i++){
    if(i%6 == 0){
      fscanf(IC, "%lf", &IC_flux[i]);
    }
    else if((i+1)%6 == 0){
      fscanf(IC, "%lf", &IC_energy[i]);
  }
    else{
      fscanf(IC, "%lf", &IC_extra[i]);
    }
  }
  */
  fclose(IC);
  fclose(brems);
  return(0);
}
