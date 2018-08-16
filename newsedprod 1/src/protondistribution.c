#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structure.h"
#include "constant.h"
#include "string.h"
#include "globalparameters.h"
#include "protondistributionlib.h"
#include "coolingprotonlib.h"
#include "pulsarmisclib.h"
double *getgamma1andGMINcontinuousproton(){
  double *gamma1GMIN;
  gamma1GMIN=malloc(2*sizeof(double));
   if(source->continuous->powerlawproton!=NULL){
     gamma1GMIN[0]=source->continuous->powerlawproton->GMIN;
     gamma1GMIN[1]=GMINproton(source->continuous->powerlawproton->GMAX);
   }

   else if(source->continuous->brokenpowerlawproton!=NULL){
     gamma1GMIN[0]=source->continuous->brokenpowerlawproton->GMIN;
     gamma1GMIN[1]=GMINproton(source->continuous->brokenpowerlawproton->GMAX);
   }

   else if(source->continuous->powerlawcutoffproton!=NULL){
    
     gamma1GMIN[0]=source->continuous->powerlawcutoffproton->GMIN;
     gamma1GMIN[1]=GMINproton(source->continuous->powerlawcutoffproton->GMAX);
   }
   
   else if(source->continuous->brokenpowerlawcutoffproton!=NULL){
     gamma1GMIN[0]=source->continuous->brokenpowerlawproton->GMIN;
     gamma1GMIN[1]=GMINproton(source->continuous->brokenpowerlawproton->GMAX);
   }

   return gamma1GMIN;
}  

double *protonGMINGMAX(){
  double *GMINMAX;
  // printf("GMAX=%.3e------\n",source->continuous->powerlawcutoffproton->GMAX);
  // sleep(3);
  GMINMAX=malloc(2*sizeof(double*));

  if(source->continuous!=NULL){
    if(source->continuous->powerlawproton!=NULL){
      if(strcmp(source->continuous->evolution,"yes")==0){
	GMINMAX[0]=GMINproton(source->continuous->powerlawproton->GMIN);
      }
      else GMINMAX[0]=source->continuous->powerlawproton->GMIN;
      GMINMAX[1]=source->continuous->powerlawproton->GMAX;
      printf("GMAX=%.3e\n",GMINMAX[1]);
    }
   else if(source->continuous->brokenpowerlawproton!=NULL){
     if(strcmp(source->continuous->evolution,"yes")==0){
	GMINMAX[0]=GMINproton(source->continuous->brokenpowerlawproton->GMIN);
      }
     else GMINMAX[0]=source->continuous->brokenpowerlawproton->GMIN;
      GMINMAX[1]=source->continuous->brokenpowerlawproton->GMAX;
     
   }
   else if(source->continuous->powerlawcutoffproton!=NULL){
     if(strcmp(source->continuous->evolution,"yes")==0){
       GMINMAX[0]=GMINproton(source->continuous->powerlawcutoffproton->GMIN);
     }
     else GMINMAX[0]=source->continuous->powerlawcutoffproton->GMIN;
      GMINMAX[1]=source->continuous->powerlawcutoffproton->GMAX;
     
    }
   else if(source->continuous->brokenpowerlawcutoffproton!=NULL){
     if(strcmp(source->continuous->evolution,"yes")==0){
       GMINMAX[0]=GMINproton(source->continuous->brokenpowerlawcutoffproton->GMIN);
      }
     else GMINMAX[0]=source->continuous->brokenpowerlawcutoffproton->GMIN;
      GMINMAX[1]=source->continuous->brokenpowerlawcutoffproton->GMAX;
      
   }
  }



  /*Now for the impulsive case*/
 if(source->impulsive!=NULL){
    if(source->impulsive->powerlawproton!=NULL){
      if(strcmp(source->impulsive->evolution,"yes")==0){
	GMINMAX[0]=GMINproton(source->impulsive->powerlawproton->GMIN);
	GMINMAX[1]=GMINproton(source->impulsive->powerlawproton->GMAX);
      }
      else {
	GMINMAX[0]=source->impulsive->powerlawproton->GMIN;
	GMINMAX[1]=source->impulsive->powerlawproton->GMAX;
      }
    }
    
    else if(source->impulsive->brokenpowerlawproton!=NULL){
     if(strcmp(source->impulsive->evolution,"yes")==0){
	GMINMAX[0]=GMINproton(source->impulsive->brokenpowerlawproton->GMIN);
	GMINMAX[1]=GMINproton(source->impulsive->brokenpowerlawproton->GMAX);
      }
     else {
       GMINMAX[0]=source->impulsive->brokenpowerlawproton->GMIN;
       GMINMAX[1]=source->impulsive->brokenpowerlawproton->GMAX;
     }
    }

    else if(source->impulsive->powerlawcutoffproton!=NULL){
     if(strcmp(source->impulsive->evolution,"yes")==0){
       
       GMINMAX[0]=GMINproton(source->impulsive->powerlawcutoffproton->GMIN);
       GMINMAX[1]=GMINproton(source->impulsive->powerlawcutoffproton->GMAX);
       printf("GMAX IS %.3e compare to %.3e\n",GMINMAX[1],source->impulsive->powerlawcutoffproton->GMAX);
     }
     else {
       GMINMAX[0]=source->impulsive->powerlawcutoffproton->GMIN;
       GMINMAX[1]=source->impulsive->powerlawcutoffproton->GMAX;
     }
     
    }
   
    else if(source->impulsive->brokenpowerlawcutoffproton!=NULL){
     if(strcmp(source->impulsive->evolution,"yes")==0){
       GMINMAX[0]=GMINproton(source->impulsive->brokenpowerlawcutoffproton->GMIN);
       GMINMAX[1]=GMINproton(source->impulsive->brokenpowerlawcutoffproton->GMAX);
      }
     else {
       GMINMAX[0]=source->impulsive->brokenpowerlawcutoffproton->GMIN;
       GMINMAX[1]=source->impulsive->brokenpowerlawcutoffproton->GMAX;
     }
   }
 }

 return GMINMAX;
}

void constantprotondistribution(){
  /* GET the protonsize*/
  int i;
  double normalisation;
  double alpha,beta,delta;
  double Ebreak,Ecutoff;
  double *GMINMAX; //pointer to the value
  int protonsize;
  GMINMAX=protonGMINGMAX();
  printf("GMIN=%.3e, GMAX=%.3e\n",GMINMAX[0],GMINMAX[1]); 
  protonsize=(int)((log10(GMINMAX[1])-log10(GMINMAX[0]))*100.0/STEP+1);  
  protonarraysize=protonsize; //pass into the globalvariable
  printf("protonarraysize=%d\n protonsize=%d\n",protonarraysize,protonsize);
  
  /*initialize the protondistribution to a double array*/
  protondistribution=malloc(2*sizeof(double*));
  protondistribution[0]=malloc(protonsize*sizeof(double));
  protondistribution[1]=malloc(protonsize*sizeof(double));
 
  /*initialisation of the protondistribution array*/
  for(i=0;i<protonsize;i++){
    protondistribution[0][i]=0;
    protondistribution[1][i]=0;
  }
 printf("protonarraysize2=%d\n",protonarraysize);
 
  /*Now we need to get the normalisation*/
  /*Case 1 continuous*/
  if(source->continuous!=NULL){
    
  normalisation=normalisationprotoncontinuous(source->continuous->time);
   printf("constantprotondistributionfunction!!!!\n");
  for(i=0;i<protonsize;i++){
    protondistribution[0][i]=GMINMAX[0]*pow(10,(i*STEP)/100.0);
  if(source->continuous->powerlawproton!=NULL){
    alpha=source->continuous->powerlawproton->inp1;
    protondistribution[1][i]=normalisation*pow(protondistribution[0][i]*protonmassTeV,-alpha);
  }
  else if(source->continuous->brokenpowerlawproton!=NULL){
    alpha=source->continuous->brokenpowerlawproton->inp1;
    beta=source->continuous->brokenpowerlawproton->inp2;
    Ebreak=source->continuous->brokenpowerlawproton->Ebreak;
    if(protondistribution[0][i]*protonmassTeV<Ebreak) protondistribution[1][i]=normalisation*pow(protondistribution[0][i]*protonmassTeV/Ebreak,-alpha);
    else protondistribution[1][i]=normalisation*pow(protondistribution[0][i]*protonmassTeV/Ebreak,-beta);
  } 
  else if(source->continuous->powerlawcutoffproton!=NULL){
  alpha=source->continuous->powerlawcutoffproton->inp1;
  beta=source->continuous->powerlawcutoffproton->expinp;
  Ecutoff=source->continuous->powerlawcutoffproton->Ecutoff;
  protondistribution[1][i]=normalisation*pow(protondistribution[0][i]*protonmassTeV,-alpha)*exp(-pow(protondistribution[0][i]*protonmassTeV/Ecutoff,beta));
  }

  else if(source->continuous->brokenpowerlawcutoffproton!=NULL){
    alpha=source->continuous->brokenpowerlawcutoffproton->inp1;
    beta=source->continuous->brokenpowerlawcutoffproton->inp2;
    delta=source->continuous->brokenpowerlawcutoffproton->expinp;
    Ebreak=source->continuous->brokenpowerlawcutoffproton->Ebreak;
    Ecutoff=source->continuous->brokenpowerlawcutoffproton->Ecutoff;
    if(protondistribution[0][i]*protonmassTeV<Ebreak) protondistribution[1][i]=normalisation*pow(protondistribution[0][i]*protonmassTeV/Ebreak,-alpha)*exp(-pow(protondistribution[0][i]*protonmassTeV/Ecutoff,delta));
    else protondistribution[1][i]=normalisation*pow(protondistribution[0][i]*protonmassTeV/Ebreak,-beta)*exp(-pow(protondistribution[0][i]*protonmassTeV/Ecutoff,delta));
  }

  }

  }

  /*Casel 2 impulsive*/
  else if (source->impulsive!=NULL){
    normalisation=normalisationprotonimpulsive();
  for(i=0;i<protonsize;i++){
    protondistribution[0][i]=GMINMAX[0]*pow(10,(i*STEP)/100.0);
  if(source->impulsive->powerlawproton!=NULL){
    alpha=source->impulsive->powerlawproton->inp1;
    protondistribution[1][i]=normalisation*pow(protondistribution[0][i]*protonmassTeV,-alpha);
  }
  else if(source->impulsive->brokenpowerlawproton!=NULL){
    alpha=source->impulsive->brokenpowerlawproton->inp1;
    beta=source->impulsive->brokenpowerlawproton->inp2;
    Ebreak=source->impulsive->brokenpowerlawproton->Ebreak;
    if(protondistribution[0][i]*protonmassTeV<Ebreak) protondistribution[1][i]=normalisation*pow(protondistribution[0][i]*protonmassTeV/Ebreak,-alpha);
    else protondistribution[1][i]=normalisation*pow(protondistribution[0][i]*protonmassTeV/Ebreak,-beta);
  } 
  else if(source->impulsive->powerlawcutoffproton!=NULL){
  alpha=source->impulsive->powerlawcutoffproton->inp1;
  beta=source->impulsive->powerlawcutoffproton->expinp;
  Ecutoff=source->impulsive->powerlawcutoffproton->Ecutoff;
  protondistribution[1][i]=normalisation*pow(protondistribution[0][i]*protonmassTeV,-alpha)*exp(-pow(protondistribution[0][i]*protonmassTeV/Ecutoff,beta));
  }

  else if(source->impulsive->brokenpowerlawcutoffproton!=NULL){
    alpha=source->impulsive->brokenpowerlawcutoffproton->inp1;
    beta=source->impulsive->brokenpowerlawcutoffproton->inp2;
    delta=source->impulsive->brokenpowerlawcutoffproton->expinp;
    Ebreak=source->impulsive->brokenpowerlawcutoffproton->Ebreak;
    Ecutoff=source->impulsive->brokenpowerlawcutoffproton->Ecutoff;
    if(protondistribution[0][i]*protonmassTeV<Ebreak) protondistribution[1][i]=normalisation*pow(protondistribution[0][i]*protonmassTeV/Ebreak,-alpha)*exp(-pow(protondistribution[0][i]*protonmassTeV/Ecutoff,delta));
    else protondistribution[1][i]=normalisation*pow(protondistribution[0][i]*protonmassTeV/Ebreak,-beta)*exp(-pow(protondistribution[0][i]*protonmassTeV/Ecutoff,delta));
  }

  }


  }

 printf("constantprotondistributionfunction!!!!\n");
  
}
  



  
 
  



void diracprotondistribution(){
  int i;
  double normalisation;
  /*This function gives the proton dirac distribution*/
  protondistribution=malloc(2*sizeof(double*));
  for(i=0;i<2;i++){
  protondistribution[i]=malloc(sizeof(double));
  }
  
  normalisation=normalisationprotondirac();
  
  protondistribution[0][0]=source->dirac->Edirac/protonmassTeV;
  protondistribution[1][0]=normalisation;
  printf("Ee dirac= %.2e Ee flux= %.3e \n",protondistribution[0][0],protondistribution[1][0]);

}




double normalisationprotoncontinuous(double time){
  double tot=0.0;
  double alpha,beta,delta;// The spectrum index
  double Ebreak,Ecutoff;
  int i;
  double gamma,dgamma;
  double normalisation;
   double age;

 if (source->continuous!=NULL) age=source->continuous->time;
 else if (source->impulsive!=NULL) age=source->impulsive->time;
   
   if (strcmp(source->continuous->spectratypeproton,"powerlaw")==0){
     alpha=source->continuous->powerlawproton->inp1;
     for(i=0;i<log10(source->continuous->powerlawproton->GMAX)*100;i++){
       gamma=pow(10,(double)i/100.0);
       dgamma=gamma*(pow(10,1.0/100.0)-1);
       tot +=pow(protonmassTeV*gamma,-alpha+1.0)*dgamma*protonmassTeV;
     }
   }

  else if (strcmp(source->continuous->spectratypeproton,"brokenpowerlaw")==0){
    alpha=source->continuous->brokenpowerlawproton->inp1;
    beta=source->continuous->brokenpowerlawproton->inp2;
    Ebreak=source->continuous->brokenpowerlawproton->Ebreak;
    for(i=0;i<log10(source->continuous->brokenpowerlawproton->GMAX)*100;i++){
       gamma=pow(10,(double)i/100.0);
       dgamma=gamma*(pow(10,1.0/100.0)-1);
       if(gamma*protonmassTeV<Ebreak) tot +=pow(protonmassTeV*gamma/Ebreak,-alpha+1)*Ebreak*dgamma*protonmassTeV;
       else tot +=pow(protonmassTeV*gamma/Ebreak,-beta+1)*Ebreak*dgamma*protonmassTeV;
     }
  }

   if (strcmp(source->continuous->spectratypeproton,"powerlawcutoff")==0){
     
     alpha=source->continuous->powerlawcutoffproton->inp1;
     beta=source->continuous->powerlawcutoffproton->expinp;
     Ecutoff=source->continuous->powerlawcutoffproton->Ecutoff;
     
       for(i=0;i<log10(source->continuous->powerlawcutoffproton->GMAX)*100;i++){
       gamma=pow(10,(double)i/100.0);
       dgamma=gamma*(pow(10,1.0/100.0)-1);
       tot +=pow(protonmassTeV*gamma,-alpha+1)*exp(-pow(protonmassTeV*gamma/Ecutoff,beta))*dgamma*protonmassTeV;
      
     }
   }


   if (strcmp(source->continuous->spectratypeproton,"brokenpowerlawcutoff")==0){
     alpha=source->continuous->brokenpowerlawcutoffproton->inp1;
     beta=source->continuous->brokenpowerlawcutoffproton->inp2;
     delta=source->continuous->brokenpowerlawcutoffproton->expinp;
     Ebreak=source->continuous->brokenpowerlawcutoffproton->Ebreak;
     Ecutoff=source->continuous->brokenpowerlawcutoffproton->Ecutoff;
     for(i=0;i<log10(source->continuous->brokenpowerlawcutoffproton->GMAX)*100;i++){
       gamma=pow(10,(double)i/100.0);
       dgamma=gamma*(pow(10,1.0/100.0)-1);
       if(gamma*protonmassTeV<Ebreak) tot +=pow(protonmassTeV*gamma/Ebreak,-alpha+1)*Ebreak*exp(-pow(protonmassTeV*gamma/Ecutoff,delta))*dgamma*protonmassTeV;
       else tot +=pow(protonmassTeV*gamma/Ebreak,-beta+1)*Ebreak*exp(-pow(protonmassTeV*gamma/Ecutoff,delta))*dgamma*protonmassTeV;
     }
   }
   /*NEED TO INPUT THE PULSAR PARAMETERS HERE*/
   if(source->continuous->pulsartype!=NULL){
     normalisation=pulsarenergy(age-time)*ergtoTeV*source->continuous->Efractionproton*(1-source->continuous->pulsartype->fractionmag)/tot;
   }
   else normalisation=source->continuous->L0*ergtoTeV*source->continuous->Efractionproton/tot;
      
   return normalisation;
}


double  pulsarenergy (double time) {
  double res;
  double P0=source->continuous->pulsartype->period*pow(1-source->continuous->pulsartype->perioddot*(source->continuous->pulsartype->brakingindex-1)*source->continuous->time/(source->continuous->pulsartype->period),1.0/(source->continuous->pulsartype->brakingindex-1));
  double Pdot0=pow(source->continuous->pulsartype->period/P0,source->continuous->pulsartype->brakingindex-2)*source->continuous->pulsartype->perioddot;
  double t0=P0/((source->continuous->pulsartype->brakingindex-1)*Pdot0);

  res = source->continuous->L0*pow(1.0+time/t0,-(source->continuous->pulsartype->brakingindex+1)/(source->continuous->pulsartype->brakingindex-1));
 
  return res;

}





double normalisationprotonimpulsive(){
double tot=0.0;
  double alpha,beta,delta;// The spectrum index
  double Ecutoff,Ebreak;
  double dgamma;
  int i;
  double gamma;
  double normalisation;
   if (strcmp(source->impulsive->spectratypeproton,"powerlaw")==0){
     alpha=source->impulsive->powerlawproton->inp1;
     for(i=0;i<log10(source->impulsive->powerlawproton->GMAX)*100;i++){
       gamma=pow(10,(double)i/100.0);
       dgamma=gamma*(pow(10,1.0/100.0)-1);
       tot +=pow(protonmassTeV*gamma,-alpha+1.0)*dgamma*protonmassTeV;
     }
   }

  else if (strcmp(source->impulsive->spectratypeproton,"brokenpowerlaw")==0){
    alpha=source->impulsive->brokenpowerlawproton->inp1;
    beta=source->impulsive->brokenpowerlawproton->inp2;
    Ebreak=source->impulsive->brokenpowerlawproton->Ebreak;
    for(i=0;i<log10(source->impulsive->brokenpowerlawproton->GMAX)*100;i++){
       gamma=pow(10,(double)i/100.0);
       dgamma=gamma*(pow(10,1.0/100.0)-1);
       if(gamma*protonmassTeV<Ebreak) tot +=pow(protonmassTeV*gamma/Ebreak,-alpha+1)*Ebreak*dgamma*protonmassTeV;
       else tot +=pow(protonmassTeV*gamma/Ebreak,-beta+1)*Ebreak*dgamma*protonmassTeV;
     }
  }

   if (strcmp(source->impulsive->spectratypeproton,"powerlawcutoff")==0){
     alpha=source->impulsive->powerlawcutoffproton->inp1;
     beta=source->impulsive->powerlawcutoffproton->expinp;
     Ecutoff=source->impulsive->powerlawcutoffproton->Ecutoff;
     for(i=0;i<log10(source->impulsive->powerlawcutoffproton->GMAX)*100;i++){
       gamma=pow(10,(double)i/100.0);
       dgamma=gamma*(pow(10,1.0/100.0)-1);
       tot +=pow(protonmassTeV*gamma,-alpha+1)*exp(-pow(protonmassTeV*gamma/Ecutoff,beta))*dgamma*protonmassTeV;
     }
   }


   if (strcmp(source->impulsive->spectratypeproton,"brokenpowerlawcutoff")==0){
     alpha=source->impulsive->brokenpowerlawcutoffproton->inp1;
     beta=source->impulsive->brokenpowerlawcutoffproton->inp2;
     delta=source->impulsive->brokenpowerlawcutoffproton->expinp;
     Ebreak=source->impulsive->brokenpowerlawcutoffproton->Ebreak;
     Ecutoff=source->impulsive->brokenpowerlawcutoffproton->Ecutoff;
     for(i=0;i<log10(source->impulsive->brokenpowerlawcutoffproton->GMAX)*100;i++){
       gamma=pow(10,(double)i/100.0);
       dgamma=gamma*(pow(10,1.0/100.0)-1);
       if(gamma*protonmassTeV<Ebreak) tot +=pow(protonmassTeV*gamma/Ebreak,-alpha+1)*Ebreak*exp(-pow(protonmassTeV*gamma/Ecutoff,delta))*dgamma*protonmassTeV;
       else tot +=pow(protonmassTeV*gamma/Ebreak,-beta+1)*Ebreak*exp(-pow(protonmassTeV*gamma/Ecutoff,delta))*dgamma*protonmassTeV;
     }
   }
   
   normalisation=source->impulsive->E0*ergtoTeV*source->impulsive->Efractionproton/tot;
     
   return normalisation;


}




double normalisationprotondirac(){
  double normalisation,total;
  /*in this function it is very brief, the distribution function resumes to a simple dirac function*/
  total=source->dirac->Edirac;
  printf("Edirac=%.3e\n",total);
  normalisation=source->dirac->E0*ergtoTeV/total;

  return normalisation;
}




void evolutionprotoncontinuous(){
  /*This function takes into account the cooling effect due to non-thermal interactions */
  int i,j;
  double *GMINMAX; //gamma2 GMAX
  int protonsize;
  double *gamma1GMIN; //gamma1 GMIN
  
  double *tmpg0reso;
  double normalisation;
  double gammastart,dgamma,gammadot;
  

  double tau;
  double dtau;
 printf("AND NOW\n"); 
  GMINMAX=protonGMINGMAX();
  printf("Gamma2=%.3e, GMAX=%.3e\n",GMINMAX[0],GMINMAX[1]);
  sleep(4);
  gamma1GMIN=getgamma1andGMINcontinuousproton();
  /*the first value would be GMIN and the second value would be gamma1*/
  protonsize=(int)((log10(GMINMAX[1])-log10(GMINMAX[0]))*100.0/STEP+1);  
  protonarraysize=protonsize;
 
  /*initialize the protondistribution to a triple array
 1/ gamma
 2/protonflux
 3/gamma_0 
 4/resolution   */
  protondistribution=malloc(4*sizeof(double*));
  for(i=0;i<4;i++){
  protondistribution[i]=malloc(protonsize*sizeof(double));
  }

  /*initialisation of the protondistribution array*/
  for(j=0;j<4;j++){
  for(i=0;i<protonsize;i++){
    protondistribution[j][i]=0;
  } 
  }
  // printf("debug1 GMIN=%.3e gamma1=%.3e\n",gamma1GMIN[0],gamma1GMIN[1]);
 
  /*Now onto the main topic, we start by the first case where gamma1>GMINMAX[0], -----MOST PROBABLE CASE-----*/
  if(gamma1GMIN[0]<gamma1GMIN[1]){
     
    for(i=0;i<protonsize;i++){
      protondistribution[0][i]=GMINMAX[0]*pow(10,(double)(i*STEP)/100.0);
      tmpg0reso=G0proton(protondistribution[0][i],source->continuous->time);
      protondistribution[2][i]=tmpg0reso[0];
      protondistribution[3][i]=tmpg0reso[1];
      free(tmpg0reso);
	printf("escape= %s\n",source->continuous->escape);
     
      /*case where the gamma and its G0 is the same*/
      if(protondistribution[2][i]==protondistribution[0][i]){
	dtau=source->continuous->time/1000.0;
	for(tau=0.0;tau<source->continuous->time;tau +=dtau){
	if(strcmp(source->continuous->escape,"yes")==0){
	  protondistribution[1][i] +=getinitialcontinuousprotondistribution(protondistribution[0][i],tau)*exp(-getescapeprotoncontinuous(protondistribution[0][i],tau))*dtau;
	  
	 
	}
	
	else protondistribution[1][i] +=getinitialcontinuousprotondistribution(protondistribution[0][i],tau)*dtau;
	  
	}
      
      }
      /*case where the resolution actually resolve g0 from gamma*/
      else{
       
	/*subcase no1 : gamma2<gamma<GMIN*/
	if((protondistribution[0][i]>=GMINMAX[0])&&(protondistribution[0][i]<gamma1GMIN[0])){
	  gammastart=gamma1GMIN[0]; //start with gamma=GMIN;
	  tau=0;  //initialise time */
	 
	  while (gammastart<protondistribution[2][i]){
	    dgamma=gammastart*(pow(10,1.0/(100*pow(10,protondistribution[3][i])))-1);
	    if(dgamma>protondistribution[2][i]-gammastart) dgamma=protondistribution[2][i]-gammastart;
	    gammadot=coolingproton(gammastart,tau);
	    tau +=dgamma/gammadot;
	    /*-----debug condition-----*/
	    if(tau>source->continuous->time) tau=source->continuous->time;
	    	if(strcmp(source->continuous->escape,"yes")==0){
		  protondistribution[1][i] +=getinitialcontinuousprotondistribution(gammastart,tau)*exp(-getescapeprotoncontinuous(gammastart,tau))*dgamma;
		}
		else {
		  protondistribution[1][i] +=getinitialcontinuousprotondistribution(gammastart,tau)*dgamma; 
		}
		gammastart +=dgamma; //gamma incrementation
	  }
	  protondistribution[1][i] /=coolingproton(protondistribution[0][i],0);
	  
	}

	/*2nd subcase: GMIN<gamma<gamma1*/
	
	if((protondistribution[0][i]>=gamma1GMIN[0])&&(protondistribution[0][i]<gamma1GMIN[1])){
	  gammastart=protondistribution[0][i]; //start with gamma;
	  tau=0;  //initialise time */
	  printf("gamma[%d]=%.3e gamma_0[%d]=%.3e GMAX=%.3e\n",i,protondistribution[0][i],i,protondistribution[2][i],GMINMAX[1]);
	  while (gammastart<protondistribution[2][i]){
	    dgamma=gammastart*(pow(10,1.0/(100*pow(10,protondistribution[3][i])))-1);
	    if(dgamma>protondistribution[2][i]-gammastart) dgamma=protondistribution[2][i]-gammastart;
	    //printf("dgamma=%.3e\n",dgamma);
	    gammadot=coolingproton(gammastart,tau);
	    tau +=dgamma/gammadot;
	    /*-----debug condition-----*/
	    if(tau>source->continuous->time) tau=source->continuous->time;
	    	if(strcmp(source->continuous->escape,"yes")==0){
		  // printf("hello? is it me you re looking for?\n");
		  protondistribution[1][i] +=getinitialcontinuousprotondistribution(gammastart,tau)*exp(-getescapeprotoncontinuous(gammastart,tau))*dgamma;
		}
		else {
		  
		  protondistribution[1][i] +=getinitialcontinuousprotondistribution(gammastart,tau)*dgamma; 
		}
		gammastart +=dgamma; //gamma incrementation
	  }
	  protondistribution[1][i] /=coolingproton(protondistribution[0][i],0);
	  
	}
	
	/* 3rd case gamma1<gamma<GMAX*/
	if((protondistribution[0][i]>=gamma1GMIN[1])&&(protondistribution[0][i]<GMINMAX[1])){
	  
	  gammastart=protondistribution[0][i]; //start with gamma;
	  tau=0;  //initialise time */
	  while (gammastart<GMINMAX[1]){
	    dgamma=gammastart*(pow(10,1.0/100.0)-1);
	    if(dgamma>GMINMAX[1]-gammastart) dgamma=GMINMAX[1]-gammastart;
	    gammadot=coolingproton(gammastart,tau);
	    tau +=dgamma/gammadot;
	    /*-----debug condition-----*/
	    if(tau>source->continuous->time) tau=source->continuous->time;
	    if(strcmp(source->continuous->escape,"yes")==0){
		  protondistribution[1][i] +=getinitialcontinuousprotondistribution(gammastart,tau)*exp(-getescapeprotoncontinuous(gammastart,tau))*dgamma;
		}
		else {
		  protondistribution[1][i] +=getinitialcontinuousprotondistribution(gammastart,tau)*dgamma; 
		}
	    	gammastart +=dgamma; //gamma incrementation
	  }
	  protondistribution[1][i] /=coolingproton(protondistribution[0][i],0);
	}
      }
    }
  }

      /*----END OF THE CASE WHERE gamma1>GMIN-------*/




      /*----2nd case : MUCH LESS LIKELY gamma1<GMIN-----*/
      else if(gamma1GMIN[1]<gamma1GMIN[0]){
	
	for(i=0;i<protonsize;i++){
      protondistribution[0][i]=GMINMAX[0]*pow(10,(double)(i*STEP)/100.0);
      tmpg0reso=G0proton(protondistribution[0][i],source->continuous->time);
      protondistribution[2][i]=tmpg0reso[0];
      protondistribution[3][i]=tmpg0reso[1];
      free(tmpg0reso);
      	
      /*case where the gamma and its G0 is the same*/
      if(protondistribution[2][i]==protondistribution[0][i]){
	dtau=source->continuous->time/1000.0;
	for(tau=0.0;tau<source->continuous->time;tau +=dtau){
	if(strcmp(source->continuous->escape,"yes")==0){
	  protondistribution[1][i] +=getinitialcontinuousprotondistribution(protondistribution[0][i],tau)*exp(-getescapeprotoncontinuous(protondistribution[0][i],tau))*dtau;
	  printf("protondistributionbis=%.3e\n",protondistribution[1][i]);
	}
	
	else protondistribution[1][i] +=getinitialcontinuousprotondistribution(protondistribution[0][i],tau)*dtau;
	  
	}
      
      }
      /*case where the resolution actually resolve g0 from gamma*/
      else{
       
	/*subcase no1 : gamma2<gamma<gamma1*/
	if((protondistribution[0][i]>=GMINMAX[0])&&(protondistribution[0][i]<gamma1GMIN[1])){
	  gammastart=gamma1GMIN[0]; //start with gamma=GMIN;
	  tau=0;  //initialise time */
	  while (gammastart<protondistribution[2][i]){
	    dgamma=gammastart*(pow(10,1.0/(100*pow(10,protondistribution[3][i])))-1);
	    if(dgamma>protondistribution[2][i]-gammastart) dgamma=protondistribution[2][i]-gammastart;
	    gammadot=coolingproton(gammastart,tau);
	    tau +=dgamma/gammadot;
	    /*-----debug condition-----*/
	    if(tau>source->continuous->time) tau=source->continuous->time;
	    	if(strcmp(source->continuous->escape,"yes")==0){
		  protondistribution[1][i] +=getinitialcontinuousprotondistribution(gammastart,tau)*exp(-getescapeprotoncontinuous(gammastart,tau))*dgamma;
		}
		else {
		  protondistribution[1][i] +=getinitialcontinuousprotondistribution(gammastart,tau)*dgamma; 
		}
		gammastart +=dgamma; //gamma incrementation
	  }
	  protondistribution[1][i] /=coolingproton(protondistribution[0][i],0);
	  
	}

	/*2nd subcase: gamma1<gamma<GMIN*/
	
	if((protondistribution[0][i]>=gamma1GMIN[1])&&(protondistribution[0][i]<gamma1GMIN[0])){
	  gammastart=gamma1GMIN[0]; //start with gamma;
	  tau=0;  //initialise time */
	  while (gammastart<GMINMAX[1]){
	    dgamma=gammastart*(pow(10,1.0/100.0)-1);
	    if(dgamma>protondistribution[2][i]-gammastart) dgamma=protondistribution[2][i]-gammastart;
	    gammadot=coolingproton(gammastart,tau);
	    tau +=dgamma/gammadot;
	    /*-----debug condition-----*/
	    if(tau>source->continuous->time) tau=source->continuous->time;
	    	if(strcmp(source->continuous->escape,"yes")==0){
		  protondistribution[1][i] +=getinitialcontinuousprotondistribution(gammastart,tau)*exp(-getescapeprotoncontinuous(gammastart,tau))*dgamma;
		}
		else {
		  protondistribution[1][i] +=getinitialcontinuousprotondistribution(gammastart,tau)*dgamma; 
		}
		gammastart +=dgamma; //gamma incrementation
	  }
	  protondistribution[1][i] /=coolingproton(protondistribution[0][i],0);
	}
	
	/* 3rd case GMIN<gamma<GMAX*/
	if((protondistribution[0][i]>=gamma1GMIN[1])&&(protondistribution[0][i]<GMINMAX[1])){
	  gammastart=protondistribution[0][i]; //start with gamma;
	  tau=0;  //initialise time */
	  while (gammastart<GMINMAX[1]){
	    dgamma=gammastart*(pow(10,1.0/100.0)-1);
	    if(dgamma>GMINMAX[1]-gammastart) dgamma=GMINMAX[1]-gammastart;
	    gammadot=coolingproton(gammastart,tau);
	    tau +=dgamma/gammadot;
	    /*-----debug condition-----*/
	    if(tau>source->continuous->time) tau=source->continuous->time;
	    if(strcmp(source->continuous->escape,"yes")==0){
		  protondistribution[1][i] +=getinitialcontinuousprotondistribution(gammastart,tau)*exp(-getescapeprotoncontinuous(gammastart,tau))*dgamma;
		}
		else {
		  protondistribution[1][i] +=getinitialcontinuousprotondistribution(gammastart,tau)*dgamma; 
		}
	    	gammastart +=dgamma; //gamma incrementation
	  }
	  protondistribution[1][i] /=coolingproton(protondistribution[0][i],0);
	}
      }
	}
      }
}




    double getinitialcontinuousprotondistribution(double gamma,double tau){
  /*this function aims at determining the initial number of proton at loretz factor gamma*/
      

  double alpha,beta,delta;
  double Ebreak,Ecutoff;
  double normalisation;
  double protondistribution2;
  
  normalisation=normalisationprotoncontinuous(tau);
  if(source->continuous->powerlawproton!=NULL){
    alpha=source->continuous->powerlawproton->inp1;
    protondistribution2=normalisation*pow(gamma*protonmassTeV,-alpha);
  }
  else if(source->continuous->brokenpowerlawproton!=NULL){
    alpha=source->continuous->brokenpowerlawproton->inp1;
    beta=source->continuous->brokenpowerlawproton->inp2;
    Ebreak=source->continuous->brokenpowerlawproton->Ebreak;
    if(gamma*protonmassTeV<Ebreak) protondistribution2=normalisation*pow(gamma*protonmassTeV/Ebreak,-alpha);
    else protondistribution2=pow(gamma*protonmassTeV/Ebreak,-beta);
  } 
  else if(source->continuous->powerlawcutoffproton!=NULL){
  alpha=source->continuous->powerlawcutoffproton->inp1;
  beta=source->continuous->powerlawcutoffproton->expinp;
  Ecutoff=source->continuous->powerlawcutoffproton->Ecutoff;
  protondistribution2=normalisation*pow(gamma*protonmassTeV,-alpha)*exp(-pow(gamma*protonmassTeV/Ecutoff,beta));
  //printf("protondistribution2=%.3e\n",protondistribution2);
  }

  else if(source->continuous->brokenpowerlawcutoffproton!=NULL){
    alpha=source->continuous->brokenpowerlawcutoffproton->inp1;
    beta=source->continuous->brokenpowerlawcutoffproton->inp2;
    delta=source->continuous->brokenpowerlawcutoffproton->expinp;
    Ebreak=source->continuous->brokenpowerlawcutoffproton->Ebreak;
    Ecutoff=source->continuous->brokenpowerlawcutoffproton->Ecutoff;
    if(gamma*protonmassTeV<Ebreak) protondistribution2=normalisation*pow(gamma*protonmassTeV/Ebreak,-alpha)*exp(-pow(gamma*protonmassTeV/Ecutoff,delta));
    else protondistribution2=normalisation*pow(gamma*protonmassTeV/Ebreak,-beta)*exp(-pow(gamma*protonmassTeV/Ecutoff,delta));
  }

  return protondistribution2;
}


double getescapeprotoncontinuous(double gamma,double tau){
  /*This function compute the escape constant and determine the amount of particles that have escaped the system*/
  double *tmpg0reso;
  double *GMINMAX;
  double escaperate;
  double Bohmcte; //bohm constant from Manolakou et al except we use a full 3D diffusion radius instead of a radial only 1D expansion of the particles, this can be changed later on */
  double totalescape=0.0;
  double gammaprime=gamma;
  double dgamma,gammadot;
  double Bfield,size;
  double tau2=0.0; //for dgamma/gammadot and get the B(t)
  double dtau2;
  if(source->continuous!=NULL){
    Bfield=source->continuous->Bfield;
    size=source->continuous->size;
  }
  else if(source->impulsive!=NULL){
    Bfield=source->impulsive->Bfield;
    size=source->impulsive->size;
  }
  Bohmcte=6*protonmassTeV*TeVtoJoules/(3*echarge*Bfield*1E-6*GausstoTesla);
  
  tmpg0reso=G0proton(gamma,tau);
  //printf("gammaprimeorigin=%.3e tmpg0res0= %.f",gammaprime,tmpg0reso[1]);
  /*case where gamma=g0*/
  if(gamma==tmpg0reso[0]){
    if (source->continuous->pulsartype==NULL) totalescape=Bohmcte/(pow(size*parsectocm,2))*gamma*tau;
    else {
      dtau2=tau/1000.0;
      for(tau2=1*yeartosec;tau2<tau;tau2 +=dtau2){

	totalescape +=Bohmcte/(pow(getsize(tau2),2))*gamma*dtau2;
      }
    //printf("totalescapeproton=%.3e\n",totalescape);
    return totalescape;
    }
  }
  else {
    if(tmpg0reso[0]>1E12) tmpg0reso[0]=1E12; //DEBUG TEST 
    /* ---NEED TO IMPOSE CONDITION ON G0>GMAX----*/
    while (gammaprime<tmpg0reso[0]){
      if (tmpg0reso[1]>=0) dgamma=gammaprime*(pow(10,1.0/(100*pow(10,tmpg0reso[1])))-1);
      else dgamma=gammaprime*(pow(10,1.0/(100.0))-1);

      if(dgamma>=tmpg0reso[0]-gammaprime) dgamma=tmpg0reso[0]-gammaprime;
      gammadot=coolingproton(gammaprime,tau2);
      tau2 +=dgamma/gammadot;
      if (source->continuous->pulsartype==NULL){
	totalescape+=Bohmcte/(pow(size*parsectocm,2))*gamma/gammadot*dgamma;
      }
      else{
	totalescape+=Bohmcte/(pow(getsize(tau2),2))*gamma/gammadot*dgamma;
      }
      gammaprime +=dgamma;
      
    }
  }
				      
  return totalescape;

}






void evolutionprotonimpulsive(){
/*This function takes into account the cooling effect due to non-thermal interactions */
  int i,j;
  double *GMINMAX; //gamma2 GMAX
  int protonsize;
  
  
  double *tmpg0reso;
  double normalisation;
  double dgamma;
  

 
  GMINMAX=protonGMINGMAX();
  printf("MOVE YOUR BODY\n");
  /*the first value would be GMIN and the second value would be gamma1*/
  protonsize=(int)((log10(GMINMAX[1])-log10(GMINMAX[0]))*100.0/STEP+1);  
  protonarraysize=protonsize;
  printf("GMIN=%.3e GMAX=%.3e protonsize=%d\n",GMINMAX[0],GMINMAX[1],protonsize);
  /*initialize the protondistribution to a quadruple array
 1/ gamma
 2/protonflux
 3/gamma_0 
 4/resolution   */
  protondistribution=malloc(4*sizeof(double*));
  for(i=0;i<4;i++){
  protondistribution[i]=malloc(protonsize*sizeof(double));
  }
  /*initialisation of the protondistribution array*/
  for(j=0;j<4;j++){
  for(i=0;i<protonsize;i++){
    protondistribution[j][i]=0;
  } 
  }

  /*In this function , there is no cases as in the continuous function because the proton distribution is dependent on a dirac function only*/
  for(i=0;i<protonsize;i++){
    protondistribution[0][i]=GMINMAX[0]*pow(10,(double)(i*STEP)/100.0);
    printf("protondistribution[0][%d]=%.3e\n",i,protondistribution[0][i]);
    tmpg0reso=G0proton(protondistribution[0][i],source->impulsive->time);
    protondistribution[2][i]=tmpg0reso[0];
    protondistribution[3][i]=tmpg0reso[1];
    free(tmpg0reso);
    printf("tmpg0reso=%.3e\n",protondistribution[2][i]);
    /*Now we got the escape option, we can use the continuous escape option that is still relevant to this particular case*/
    if(strcmp(source->impulsive->escape,"yes")==0){
      printf("protonimpulsive?\n");
      printf("time=%.3e\n",source->impulsive->time);
      protondistribution[1][i] =getinitialimpulsiveprotondistribution(protondistribution[2][i])*exp(-getescapeprotoncontinuous(protondistribution[0][i],source->impulsive->time))*dgamma;
      
    }
    else {
      protondistribution[1][i] =getinitialimpulsiveprotondistribution(protondistribution[2][i]); 
      // printf("protonimpulsive?\n");
    }
    protondistribution[1][i] *=coolingproton(protondistribution[2][i],source->impulsive->time)/coolingproton(protondistribution[0][i],0);
     printf("protondistribution[1][%d]=%.3e\n",i,protondistribution[1][i]);
  }

}								       
















 double getinitialimpulsiveprotondistribution(double gamma){
  /*this function aims at determining the initial number of proton at loretz factor gamma*/


  double alpha,beta,delta;
  double Ebreak,Ecutoff;
  double normalisation;
  double protondistribution2;
  normalisation=normalisationprotonimpulsive();
  if(source->impulsive->powerlawproton!=NULL){
    alpha=source->impulsive->powerlawproton->inp1;
    protondistribution2=normalisation*pow(gamma*protonmassTeV,-alpha);
  }
  else if(source->impulsive->brokenpowerlawproton!=NULL){
    alpha=source->impulsive->brokenpowerlawproton->inp1;
    beta=source->impulsive->brokenpowerlawproton->inp2;
    Ebreak=source->impulsive->brokenpowerlawproton->Ebreak;
    if(gamma*protonmassTeV<Ebreak) protondistribution2=normalisation*pow(gamma*protonmassTeV/Ebreak,-alpha);
    else protondistribution2=pow(gamma*protonmassTeV/Ebreak,-beta);
  } 
  else if(source->impulsive->powerlawcutoffproton!=NULL){
  alpha=source->impulsive->powerlawcutoffproton->inp1;
  beta=source->impulsive->powerlawcutoffproton->expinp;
  Ecutoff=source->impulsive->powerlawcutoffproton->Ecutoff;
  protondistribution2=normalisation*pow(gamma*protonmassTeV,-alpha)*exp(-pow(gamma*protonmassTeV/Ecutoff,beta));
  printf("are we here?\n");
  }

  else if(source->impulsive->brokenpowerlawcutoffproton!=NULL){
    alpha=source->impulsive->brokenpowerlawcutoffproton->inp1;
    beta=source->impulsive->brokenpowerlawcutoffproton->inp2;
    delta=source->impulsive->brokenpowerlawcutoffproton->expinp;
    Ebreak=source->impulsive->brokenpowerlawcutoffproton->Ebreak;
    Ecutoff=source->impulsive->brokenpowerlawcutoffproton->Ecutoff;
    if(gamma*protonmassTeV<Ebreak) protondistribution2=normalisation*pow(gamma*protonmassTeV/Ebreak,-alpha)*exp(-pow(gamma*protonmassTeV/Ecutoff,delta));
    else protondistribution2=normalisation*pow(gamma*protonmassTeV/Ebreak,-beta)*exp(-pow(gamma*protonmassTeV/Ecutoff,delta));
  }
  printf("protondistribution2=%.3e\n",protondistribution2);
  return protondistribution2;
 }
