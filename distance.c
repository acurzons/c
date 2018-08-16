/*
 * distance.c
 *
 * Outputs distance in deg between two RA, Dec positions 
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159

int main(int argc, char **argv)
{
  double ra1, dec1, ra2, dec2, ra0, dec0;
  double phi, t1, z1, t2, z2, az1, az2, d;

  if (argc < 4) {
    printf("Usage: Non Noff alpha\n");
    exit(0);
  }

  ra1  = atof(argv[1]);
  dec1 = atof(argv[2]);
  ra2  = atof(argv[3]);
  dec2 = atof(argv[4]); 

  ra0  = atof(argv[5]);
  dec0 = atof(argv[6]);

  phi = (dec0 - 90.0) * PI/180.0; /* Lat of observation here is set to force
				     the original source at elevation 0.0 deg.*/
  
  t1 = (ra0 - ra1) * 15.0 * PI/180.0; /*Set hours angle re: Original RA*/
  z1 = acos(sin(phi)*sin(dec1*PI/180.0) + cos(phi)*cos(dec1*PI/180.0)*
	    cos(t1));
  az1 = 180.0 * asin(cos(dec1*PI/180.0)*sin(t1)/sin(z1))/PI;
  z1 = 180.0 * z1 / PI;

  t2 = (ra0 - ra2) * 15.0 * PI/180.0;  /*Set hours angle re: Original RA*/
  z2 = acos(sin(phi)*sin(dec2*PI/180.0) + cos(phi)*cos(dec2*PI/180.0)*
	    cos(t2));
  az2 = 180.0 * asin(cos(dec2*PI/180.0)*sin(t2)/sin(z2))/PI;
  z2 = 180.0 * z2 / PI;

  d = sqrt( (z1-z2)*(z1-z2) + (az1-az2)*(az1-az2) );


 printf("%f\n", d);

}






