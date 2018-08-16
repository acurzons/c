#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define TeVtoJoules  1.60217657E-7 // eV to J conversion
#define eVtoJoules 1.60217657E-19 // eV to J conversion
#define electronmassTeV 0.510998902E-6 // in TeV
#define UCMBTeV 0.256E-12 //in TeV/cm^-3
#define PLANCK 6.62606957E-34
#define PI 3.14159265359
#define lightspeedcm 2.99792458E10 
#define thompsonXS 6.6524587E-25 // in cm^2
#define BOLTZMANN 1.380658E-23 //boltzmann factor in J/K
/*----------------------------------------------------------------------------
TESTIC.c: this program aims at testing the different regime from the inverse compton

-----------------------------------------------------------------------------*/
double inversecomptonfunction(double Ee, double Egamma, double initphotonenergy);
double blackbodydistribution2(double Urad, double Egamma, double T);
int main(int argc,char **argv){
  /*The program consists in injecting a monochromatic distribution of electron, the program should list*/
  FILE *ofile;
  char filename[100];
  
  double inversecompton;
  double temperature;
  double Energyinitphoton;
  double dEnergyinit;
  double electronenergy;
  double Egamma;//TeV
  int i,j;
  const double electronnumber=1E46;//Number of electron

  temperature=atof(argv[2]);
  electronenergy=atof(argv[3]);
  Egamma=atof(argv[4]);
  /*We need to loop over the different gamma*/
 

  
    sprintf(filename,"%s.txt",argv[1]);
    if((ofile=fopen(filename,"w+"))==NULL){
      printf("unable to write file\n");
      exit(1);
    }
  
      Energyinitphoton=2.8*BOLTZMANN*temperature/TeVtoJoules/1000.0;
      dEnergyinit=2.8*BOLTZMANN*temperature/TeVtoJoules/100.0;
      while (Energyinitphoton<2.8*BOLTZMANN*temperature/TeVtoJoules*10){
	// inversecompton=3.0/4.0*thompsonXS*lightspeedcm*pow(electronmassTeV/electronenergy,2)*blackbodydistribution2(UCMBTeV*1E12,targetenergy[0],temperature)/targetenergy[0]*inversecomptonfunction(electronenergy,Egamma,targetenergy[0])*electronnumber;
	 inversecompton=inversecomptonfunction(electronenergy,Egamma,Energyinitphoton);
	 
      
 
	 fprintf(ofile,"%.3e %.3e \n",Energyinitphoton,inversecompton);
	 Energyinitphoton +=dEnergyinit;
      }
    
      fclose(ofile);
    

  return 0;
  }

double blackbodydistribution2(double Urad, double Egamma, double T){
  double distribution2;
  Urad *=eVtoJoules;
  distribution2=15*Urad/(pow(PI*BOLTZMANN*T,4))*pow(Egamma*TeVtoJoules,2)*1.0/(exp(Egamma*TeVtoJoules/(BOLTZMANN*T))-1.0);
  distribution2 *=TeVtoJoules; //in ph cm-3 TeV-1 

  return distribution2;
}


double inversecomptonfunction(double Ee, double Egamma, double initphotonenergy){
  double inversecompton;
  double p;
  double q; //see Meyer et al for variable name ->NON DIMENSIONAL 

  // q=Egamma/(4*initphotonenergy*pow(Ee/electronmassTeV,2.0)*(1-Egamma/Ee));
  p = 4.0*initphotonenergy*Ee/pow(electronmassTeV, 2.0);
  q = Egamma/(p*(Ee-Egamma));


  
  //inversecompton=2.0*q*log(q)+(1.0+2.0*q)*(1.0-q)+1.0/2.0*(pow(4.0*initphotonenergy*Ee/pow(electronmassTeV,2)*q,2)/(1.0+4.0*initphotonenergy*Ee/pow(electronmassTeV,2)*q)*(1.0-q));

  inversecompton=2.0*q*log(q) + (1.0+2.0*q)*(1.0-q)+ 0.5*(p*q*p*q)*(1.0-q)/(1.0+p*q);
  /*if (Egamma>20 &&Egamma<30) {
    printf("inverse compton=%.3f q=%.3f p =%.3f\n",inversecompton,q,p);
    sleep(1);
    }*/
  if ((q>1)&&(q<0)){

    // printf("q=%.3e Egamma=%.3e Ee=%.3e initphoton=%.3e\n",q,Egamma,Ee,initphotonenergy);
    inversecompton=0.0;
  }
  if (inversecompton<0) inversecompton=0;
  //if(q<0) printf("IMPOSSIIIIIIIIIBLEEEEE\n");
  return inversecompton;

}
