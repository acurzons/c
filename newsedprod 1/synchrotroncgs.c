#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_sf_synchrotron.h>


#define lightspeedcm 2.99792458E10 
#define Bfield 3.0 //microGauss
#define electronmassTeV 0.510998902E-6 // in TeV
#define ergtoTeV 0.624150934
#define totalenergy 1E48 //erg
#define PLANCKcgs  6.6260755E-27 //erg s
#define parsectocm 3.0856E+18     /* cm */
#define echargecgs 4.8066E-10     /* statcoul  */
#define electronenergy 1.0
#define PI 3.14159265359
#define distance 1000.0 //pc  
#define thompsonXS 6.6524587E-25 // in cm^2
int main(){
  double totalsynchrotron=0.0;
  double Egamma;
  double dEgamma;
  double synchrotronflux;
  double synchrotronemissivity;
  int i;
  double x;
  double frequency, criticalfrequency;
  double Bessel;
  FILE *ofile;
  double totalloss;

  if ((ofile=fopen("syncanalytic.txt","w+"))==NULL){
    printf("unable to open file \n");
    exit(1);
  }

  for (i=-2000;i<800;i++){
    Egamma=pow(10,(double)i/100.0);
    dEgamma=Egamma*(pow(10,1.0/100.0)-1);
 
    /*get the frequency and critical frequency*/
    frequency=Egamma/(ergtoTeV*PLANCKcgs);
    criticalfrequency=3.0/(4.0*PI)*sqrt(2.0/3.0)*pow(electronenergy/electronmassTeV,2.0)*echargecgs*Bfield*1E-6*pow(lightspeedcm,1.0)/(electronmassTeV/ergtoTeV);
    x=frequency/criticalfrequency;
    if ((x >= 0.0) && (x<-8.0*GSL_LOG_DBL_MIN/7.0)) { 
    Bessel=gsl_sf_synchrotron_1(x);
    
    if (Bessel<1E-9) Bessel=0.0;
  }
    else Bessel=0.0;
    /*the 3.0/2.0 at the beginning represents the sin^(alpha) from the ratio emitted to received (see Blumenthal and Gould  1970)*/
    synchrotronflux=sqrt(3)*pow(echargecgs,3.0)*Bfield*1E-6*sqrt(2.0/3.0)/(electronmassTeV*Egamma*PLANCKcgs/pow(ergtoTeV,2))*Bessel/ergtoTeV*totalenergy*ergtoTeV/electronenergy; // ph TeV -1
    synchrotronemissivity=synchrotronflux/(4*PI);
    totalsynchrotron +=Egamma*synchrotronflux*dEgamma;
    fprintf(ofile,"%.3e %.3e %.3e %.3e\n",Egamma,synchrotronflux, Egamma*Egamma*1.0/ergtoTeV*synchrotronflux/(4*PI*pow(distance*parsectocm,2.0)),Egamma*Egamma*1.0/ergtoTeV*synchrotronemissivity/(4*PI*pow(distance*parsectocm,2.0)));
  }
  totalloss=4.0/3.0*thompsonXS*pow(electronenergy/electronmassTeV,2.0)*lightspeedcm*pow(Bfield*1E-6,2.0)/(8*PI)*totalenergy*ergtoTeV/electronenergy;
  printf("synchrotrontotal=%.3e totalloss=%.3e \n",totalsynchrotron,totalloss);
  fclose(ofile);
return 0;

}
