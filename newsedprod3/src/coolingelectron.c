#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structure.h"
#include "constant.h"
#include "globalparameters.h"
#include <string.h>
#include "coolingelectronlib.h"
#include "pulsarmisclib.h"



double GMINelectron(double gamma){

  double resolution=0.0;
  double startgamma, dgamma; //variable for integration purposes
  double tau=0.0, currentime=0; //variable for integration
  double gammadot; //variable for cooling
  double time;
  double finalgamma;
  printf("gammatoreduce=%.3e\n",gamma);
  if (source->continuous!=NULL){
    time=source->continuous->time;
  }
  else if (source->impulsive!=NULL){
    time=source->impulsive->time;
  }
  printf("time=%.3e\n",time);
  while((tau+currentime<time)||(((tau+currentime)/time -1>0.0001)&&(resolution<6))){ /* conditions to accept the output as relevant*/
    resolution++;
    /*initialization time variable*/
    currentime=0;
    tau=0.0;
  startgamma=gamma;
  while (tau+currentime<time){
    dgamma=startgamma*(1-pow(10,-1.0/(100.0*pow(10,resolution))));
    startgamma-=dgamma;
    currentime +=tau;
    gammadot=coolingelectron(startgamma,currentime);
    tau=dgamma/gammadot;//computation of the requested time;
   

    if (startgamma<1E-2) {
      finalgamma=1.0;
      return finalgamma;
    }
 //printf("tau=%.3e\n dgamma=%.3e gammadot=%.3e time=%.3e",currentime,dgamma,gammadot,time);
    
  }
  printf("resogamma=%.3e\n",startgamma);
  }
  

  if(resolution>=6)finalgamma=gamma;
  else finalgamma=startgamma;





  return finalgamma;
}


double *G0electron(double gamma,double time){

  double resolution=0.0;
  double startgamma, dgamma; //variable for integration purposes
  double tau=0.0, currentime=0; //variable for integration
  double gammadot; //variable for cooling
  double *finalg0reso;
  int count=0;
  
  finalg0reso=malloc(2*sizeof(double));
 

  //while((tau+currentime<time)||(((tau+currentime)/time -1>0.001)&&(resolution<5))){ /* conditions to accept the output as relevant*/
  /* I CHANGED THE COUNT FROM 3000 to 5000 and resolution 5 to 4*/
  while((tau+currentime<time)||((count<3000)&&(resolution<4))){
    resolution++;
    count=0;
    /*initialization time variable*/
    currentime=0;
    tau=0.0;
  startgamma=gamma;
  dgamma=0.0;
  while (tau+currentime<time){
    startgamma+=dgamma;
    dgamma=startgamma*(pow(10,1.0/(100.0*pow(10,resolution)))-1);
    
    currentime +=tau;
    gammadot=coolingelectron(startgamma,currentime);
    tau=dgamma/gammadot; //computation of the requested time;
   
    count++;
    if (startgamma>1E50){ //special condition in case we get past the critical value
      finalg0reso[0]=-1;
      finalg0reso[1]=0;
      return finalg0reso;
    }
  }

  }
  

  if(resolution>=4){
    finalg0reso[0]=gamma;
    finalg0reso[1]=resolution;
  }
    else {
      finalg0reso[0]=startgamma;
      finalg0reso[1]=resolution;
  
    }


  return finalg0reso;
}








double coolingelectron(double gamma,double tau) {

  /* First we declare all the constant from Manolakou 2006*/
  //check if everything is concording in unit
  double bs, bic, bc, bb, b0c, b0b, bpp,xs;
  double res;
  double Bfield,density,TempIR,energyIR;
  double bIR; //for IR sources
  double age;
  
   if (source->continuous!=NULL){
     Bfield=source->continuous->Bfield;
     density=source->continuous->density;
     TempIR=source->continuous->TempIR;
     energyIR=source->continuous->EnergyIR;
     age=source->continuous->time;
     if(source->continuous->pulsartype==NULL) bs= 1.292E-15*pow(Bfield*1.0E-3,2);
   else bs=1.292E-15*pow(getmagneticfield(age)*1.0E-3,2);
   }
  else if (source->impulsive!=NULL){
    Bfield=source->impulsive->Bfield;
    density=source->impulsive->density;
    TempIR=source->impulsive->TempIR;
    energyIR=source->impulsive->EnergyIR;
    age=source->impulsive->time;
    bs= 1.292E-15*pow(Bfield*1.0E-3,2);
  }
   
  bic=5.204E-20*(UCMBTeV*1E12);
  bIR=5.204E-20*(energyIR);
  bc= 1.491E-14*density;
  bb= 1.37E-16*density;
  b0c=-log(density)+73.4;
  b0b=log(2)-1.0/3.0;
 
  // bpp=5.98E-16*np;
   
      res= bs*pow(gamma,2)+bic*pow(gamma,2)*kleinnishimaelectron(gamma,TEMPCMB)+bIR*pow(gamma,2)*kleinnishimaelectron(gamma,TempIR)+bc*(log(gamma)+b0c)+bb*gamma*(log(gamma)+b0b) ;
   
  return res;
}






double kleinnishimaelectron( double  gamma,double temperature ) {

// From Troy Porter PhD thesis Eq 5.24...
  // Isotropic soft target photons are assumed
  // expect energies in TeV units
  double F_kn;
  
  
  F_kn=1.0*pow((1+(4*gamma*2.8*BOLTZMANN*temperature)/(electronmassTeV*TeVtoJoules)),-1.5);
  
  return F_kn;
}







