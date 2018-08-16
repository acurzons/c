#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globalparameters.h"
#include "structure.h"
#include "constant.h"
#include "coolingprotonlib.h"
#include <string.h>
#include "pulsarmisclib.h"



double GMINproton(double gamma){

  double resolution=0.0;
  double startgamma, dgamma; //variable for integration purposes
  double tau=0.0, currenttime=0; //variable for integration
  double gammadot; //variable for cooling
  double time;
  double finalgamma;
  
  if (source->continuous!=NULL){
    time=source->continuous->time;
  }
  else if (source->impulsive!=NULL){
    time=source->impulsive->time;
  }
  
  while((tau+currenttime<time)||((tau+currenttime)/time -1>0.0001)||resolution<6){ /* conditions to accept the output as relevant*/
    resolution++;
    /*initialization time variable*/
    currenttime=0;
    tau=0.0;
  startgamma=gamma;
 
  while (tau+currenttime<time){
    dgamma=startgamma*(1-pow(10,-1.0/(100.0*pow(10,resolution))));
    startgamma-=dgamma;
    currenttime +=tau;
    gammadot=coolingproton(startgamma,currenttime);
     
    // printf("cooling %.3e\n",gammadot);
    tau=dgamma/gammadot; //computation of the requested time;
    //printf("the workd\n");
    if (startgamma<1E-2){
      finalgamma=1.0;
      return finalgamma;
	}
  }

  }
  

  if(resolution>=5)finalgamma=gamma;
  else finalgamma=startgamma;




 
  return finalgamma;
}


double *G0proton(double gamma,double time){

  double resolution=0.0;
  double startgamma, dgamma; //variable for integration purposes
  double tau=0.0, currenttime=0; //variable for integration
  double gammadot; //variable for cooling
  double *finalg0reso;
  
  
  finalg0reso=malloc(2*sizeof(double));
 

  while((tau+currenttime<time)||(((tau+currenttime)/time -1>0.01)&&(resolution<5))){ /* conditions to accept the output as relevant*/
    resolution++;
    /*initialization time variable*/
    currenttime=0;
    tau=0.0;
  startgamma=gamma;
  dgamma=0.0;
  while (tau+currenttime<time){
    startgamma+=dgamma;
    dgamma=startgamma*(pow(10,1.0/(100.0*pow(10,resolution)))-1);
    
    currenttime +=tau;
    gammadot=coolingproton(startgamma,currenttime);
    tau=dgamma/gammadot; //computation of the requested time;
   

    if (startgamma>1E50){ //special condition in case we get past the critical value
      finalg0reso[0]=-1;
      finalg0reso[1]=0;
      return finalg0reso;
    }
  }

  }
  

  if(resolution>=5){
    finalg0reso[0]=gamma;
    finalg0reso[1]=resolution;
  }
    else {
      finalg0reso[0]=startgamma;
      finalg0reso[1]=resolution;
    }




  return finalg0reso;
}








double coolingproton(double gamma,double tau) {

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
     if(source->continuous->pulsartype==NULL) bs= 2.25E-22*pow(Bfield*1.0E-3,2);
     else bs=2.25E-22*pow(getmagneticfield(age)*1.0E-3,2);
   }
  else if (source->impulsive!=NULL){
    
    Bfield=source->impulsive->Bfield;
    density=source->impulsive->density;
    TempIR=source->impulsive->TempIR;
    
    energyIR=source->impulsive->EnergyIR;
    age=source->impulsive->time;
    bs= 2.25E-22*pow(Bfield*1.0E-3,2);
    
  }
   
   // printf("are we here?\n");
   
   
  bic=9.050E-30*(UCMBTeV*1E12);
  
  bIR=9.050E-30*(energyIR);
  bc= 4.27E-21*density;
  bb= 3.93E-23*density;
  b0c=-log(density)+73.4;
  b0b=log(2)-1.0/3.0;
  xs=(34.3+1.88*log(gamma*9.38E-4)+0.25*pow(log(gamma*9.38E-4),2))*1E-27;
  bpp=0.5*3E10*xs*density;
 
    if (gamma*9.38E-4<=0.00121){
 
      res= bs*pow(gamma,2)+bic*pow(gamma,2)*kleinnishimaproton(gamma,TEMPCMB)+bIR*pow(gamma,2)*kleinnishimaproton(gamma,TempIR)+bc*(log(gamma)+b0c)+bb*gamma*(log(gamma)+b0b) ;
      
	
    }
    else {
      res= bs*pow(gamma,2)+bic*pow(gamma,2)*kleinnishimaproton(gamma,TEMPCMB)+bIR*pow(gamma,2)*kleinnishimaproton(gamma,TempIR)+bc*(log(gamma)+b0c)+bb*gamma*(log(gamma)+b0b)+bpp*gamma ;
    }
 
    
  
  return res;
}






double kleinnishimaproton( double  gamma,double temperature ) {

// From Troy Porter PhD thesis Eq 5.24...
  // Isotropic soft target photons are assumed
  // expect energies in TeV units
  double F_kn;
  
  
  F_kn=1.0*pow((1+(4*gamma*2.8*BOLTZMANN*temperature)/(protonmassTeV*TeVtoJoules)),-1.5);
  
  return F_kn;
}






