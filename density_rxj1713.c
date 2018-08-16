#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"


int main(){

  FILE *rxj1713_density;

  double column_density[29];
  double density[29];
  double pctocm = 3.085677581e18; //pc to cm
  double depth[4] = {4,5,6,7}; //how many regions deep
  double region_depth = 3.2; //depth per region in pc
  double sum_density = 0.0;
  double average_density;
  int i;

  //Regional column density in units of cm^-2/1e21
  column_density[0] = 7.23429;
  column_density[1] = 7.37471;
  column_density[2] = 9.8519;
  column_density[3] = 7.77038;
  column_density[4] = 4.84366;
  column_density[5] = 6.69706;
  column_density[6] = 4.90135;
  column_density[7] = 7.04167;
  column_density[8] = 10.5121;
  column_density[9] = 10.7982;
  column_density[10] = 7.03478;
  column_density[11] = 5.94348;
  column_density[12] = 4.24037;
  column_density[13] = 5.24132;
  column_density[14] = 7.5758;
  column_density[15] = 10.4733;
  column_density[16] = 5.44116;
  column_density[17] = 4.71859;
  column_density[18] = 4.51636;
  column_density[19] = 5.09396;
  column_density[20] = 10.9893;
  column_density[21] = 12.4869;
  column_density[22] = 5.43457;
  column_density[23] = 4.94197;
  column_density[24] = 5.41388;
  column_density[25] = 5.9631;
  column_density[26] = 11.2197;
  column_density[27] = 6.7807;
  column_density[28] = 9.15628;


  rxj1713_density = fopen("rxj1713_regional_densities.txt", "w");

  //Convert regional depth to cm depth
  for(i=0; i<4; i++) depth[i] *= region_depth * pctocm;
 
  for(i=0; i<29; i++){
    //Convert column density to units of cm^-2
    column_density[i] *= 1e21;
    
    //Calculate the volume density
    if(i==0 || i==3 || i==4||i==9||i==21||i==22||i==27||i==28) density[i] = column_density[i]/depth[0];
    if(i==1 || i==2 || i==5||i==10||i==15||i==16||i==23) density[i] = column_density[i]/depth[1];
    if(i==12 || i==13)  density[i] = column_density[i]/depth[2];
    else density[i] = column_density[i]/depth[3];

    //find the average density
    sum_density += density[i];
    average_density = sum_density/29.0;

    fprintf(rxj1713_density, "%d %.3e %.3e\n", i+1, column_density[i], density[i]);
  }

  printf("Average density = %.3e \n", average_density);


  fclose(rxj1713_density);
  return(0);
}
