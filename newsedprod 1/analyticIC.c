#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TEMPCMB 2.72 // in K
#define parsectocm 3.1E18 //cm
#define BOLTZMANN 1.380658E-23 //boltzmann factor in J/K
#define eVtoJoules 1.60217657E-19 // eV to J conversion
#define TeVtoJoules  1.60217657E-7 // eV to J conversion
#define lightspeedcm 2.99792458E10 
#define electronmassTeV 0.510998902E-6 // in TeV
#define thompsonXS 6.6524587E-25 // in cm^2
#define UCMBTeV 0.256E-12 //in TeV/cm^-3
#define PI 3.14159265359
#define ergtoTeV 0.624150934
int main(int argc, char **argv){
  /*This program will compute the analytical solution for IC at low energy (thompson regime) for a planck distribution of p*/ 
  double ICplanck;
  double ICdirac;
  FILE *ofile;
  double Energyelectron;
  double totalenergy;
  int i;
  double Egamma;
  double distance;

  totalenergy=atof(argv[1]);
  Energyelectron=atof(argv[2]);
  distance=atof(argv[3]);
  if ((ofile=fopen("analytic.txt","w+"))==NULL){
    printf("Unable to write file\n");
    exit(1);
  }
  
  for(i=-1300;i<500;i++){
    Egamma=pow(10,(double)i/100.0);
    ICplanck=3.0/4.0*thompsonXS*pow(electronmassTeV/Energyelectron,2.0)*lightspeedcm*15.0/6.0*UCMBTeV/(pow(PI*BOLTZMANN*TEMPCMB/TeVtoJoules,2.0))*totalenergy*ergtoTeV/Energyelectron*pow(Egamma,2.0)/(4*PI*pow(distance*parsectocm,2.0));

    ICdirac=3.0/4.0*thompsonXS*pow(electronmassTeV/Energyelectron,2.0)*lightspeedcm*UCMBTeV/(pow(2.8*BOLTZMANN*TEMPCMB/TeVtoJoules,2.0))*totalenergy*ergtoTeV/Energyelectron*pow(Egamma,2.0)/(4*PI*pow(distance*parsectocm,2.0));


    fprintf(ofile, "%.3e %.3e %.3e \n", Egamma,ICplanck/ergtoTeV,ICdirac/ergtoTeV);
  
  }

  fclose(ofile);
  

  return 0;
}
