#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"


int main(){

  FILE *ATCA_flux;

  double flux_J[29];
  double flux_cgs[29];
  double nu = 1.357e9;
  double Energy;

  int pixelsperbeam = 92;
  int i;

  //Regional flux taken from fits file
  flux_J[0] = 15.6342;
  flux_J[1] = 25.8488;
  flux_J[2] = 26.3648;
  flux_J[3] = 64.0503;
  flux_J[4] = 2.22941;
  flux_J[5] = 23.545;
  flux_J[6] = 14.9111;
  flux_J[7] = 15.3448;
  flux_J[8] = 58.0416;
  flux_J[9] = 124.541;
  flux_J[10] = 14.6148;
  flux_J[11] = 19.2673;
  flux_J[12] = 38.5926;
  flux_J[13] = 15.5326;
  flux_J[14] = 41.1254;
  flux_J[15] = 100.752;
  flux_J[16] = 18.1381;
  flux_J[17] = 21.7754;
  flux_J[18] = 34.7982;
  flux_J[19] = 47.8357;
  flux_J[20] = 30.8126;
  flux_J[21] = 47.6862;
  flux_J[22] = 7.72616;
  flux_J[23] = 17.0534;
  flux_J[24] = 39.592;
  flux_J[25] = 18.7065;
  flux_J[26] = 23.8323;
  flux_J[27] = 5.3525;
  flux_J[28] = 0.96626;

  ATCA_flux = fopen("../newsedprod/OUTPUT/authorsplots/ATCA_rxj1713regionalflux.txt", "w");

  //Energy in TeV (E=h*nu in units of joules and then convert to eV and TeV)
  Energy = 6.626e-34 * nu / 1.6e-19 / 1e12;
 
  for(i=0; i<29; i++){
    //Convert flux in Jy/beam to erg/cm^2/s
    flux_cgs[i] = flux_J[i] * nu * 1e-23 / pixelsperbeam;
    fprintf(ATCA_flux, "%d %.3e %.3e %.3e\n", i+1, Energy, flux_J[i], flux_cgs[i]);
  }

  fclose(ATCA_flux);
  return(0);
}
