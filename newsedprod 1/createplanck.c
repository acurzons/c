#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define lightspeedcm 2.99792458E10 
#define lightspeedSI 2.99792458E8
#define protonmassTeV 938.272E-6 // in TeV
#define electronmassTeV 0.510998902E-6 // in TeV
#define pionmassTeV 134.9766E-6
#define chargedpionmassTeV 134.57018E-6
#define neutrinomassTeV 0.320E-12
#define eVtoJoules 1.60217657E-19 // eV to J conversion
#define TeVtoJoules  1.60217657E-7 // eV to J conversion

#define BOLTZMANN 1.380658E-23 //boltzmann factor in J/K
#define TEMPCMB 2.728 // in K
#define UCMBTeV 0.262E-12 //in TeV/cm^-3
#define STEP 5.0  //step between the proton/electron distribution

#define ergtoTeV 0.624150934
#define echarge 1.6E-19 // electric charge in J
#define GausstoTesla 1E-4;



#define PLANCK 6.62606957E-34
#define PI 3.14159265359
#define vacuumpermittivity  8.854187817E-12
double blackbodydistribution(float T,float Egamma);
double blackbodydistribution2(double Urad,double T,double Egamma);
int main(int argc, char **argv){
  float planckdist;
  int i;
  FILE *ofile;
  FILE *file2;
  char filename[50];
  float photonenergy;
  float temp;
  float planckdist2;
  float total,dE;
  float dirac;
  total=0.0;
  temp=atof(argv[2]);
   sprintf(filename,"%s",argv[1]);
  if((ofile=fopen(filename,"w+"))==NULL){
    printf("unable to open file \n");
    exit(1);
  }
 if((file2=fopen("logplanck.txt","w+"))==NULL){
    printf("unable to open file \n");
    exit(1);
  }

  for(i=-200;i<80;i++){
    photonenergy=pow(10,i/10.0);
    dE=photonenergy*(pow(10,1.0/10.0)-1);
    planckdist=blackbodydistribution(temp,photonenergy);
    planckdist2=blackbodydistribution2(UCMBTeV*1E12,photonenergy,temp);
    fprintf(ofile,"%.3e %.3e %.3e\n",photonenergy,planckdist,planckdist2);
    //total +=photonenergy*planckdist2*1E12*dE;
    total +=planckdist2*photonenergy*1E12*dE;
    fprintf(file2,"%d %.3e \n",i,total);
   
    
  }
 
  for (i=1;i<100;i++){
    dirac=blackbodydistribution2(UCMBTeV*1E12,(float)(i/10.0)*BOLTZMANN*temp/TeVtoJoules,temp)*(float)(i/10.0)*BOLTZMANN*temp/TeVtoJoules*(pow(10,1.0/10.0)-1);
 
    printf("total= %.3e %.3e %.2f  ph /cm -3 \n",total, dirac,i/10.0);
  }
   fclose(ofile);
   fclose(file2);
  return 0;
}


double blackbodydistribution(float T, float Egamma){
  //returns number of the planckian distribution function describing the number of photons between Egamma and Egamma+dEgamma
  double logdistribution;
  double distribution;
  

  logdistribution=(double)(log(8.0)+log(PI)+2*log(Egamma*TeVtoJoules)-3*log(PLANCK)-3*log(lightspeedcm)-log(exp(Egamma*TeVtoJoules/(BOLTZMANN*T))-1.0)+log(TeVtoJoules));
  distribution=exp(logdistribution);
  
  // distribution /=TeVtoJoules; //in ph cm-3 TeV-1 /1E30

  return distribution;
  }


double blackbodydistribution2(double Urad, double Egamma, double T){
  double distribution2;
  Urad *=eVtoJoules;
  distribution2=15*Urad/(pow(PI*BOLTZMANN*T,4))*pow(Egamma*TeVtoJoules,2)*1.0/(exp(Egamma*TeVtoJoules/(BOLTZMANN*T))-1.0);
  distribution2 *=TeVtoJoules; //in ph cm-3 TeV-1 

  return distribution2;
}
