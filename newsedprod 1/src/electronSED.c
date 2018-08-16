#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_sf_synchrotron.h>
#include "constant.h"
#include "structure.h"
#include "globalparameters.h"
#include "electronSEDlib.h"
#include "bremss.h"
#include "pulsarmisclib.h"

void electronSEDprimary(){
  int i,j,l;
  double dEgamma;
  double dEe,Ee;
  
  double EGMIN;
  
  int sizegamma;
  
  double Bfield;
  double TempIR;
  double density;
  double dEnergyinit; //in TeV
  double EnergyinitphotonCMB,Energyinitphoton;
  double IC1,ICR,ICO;
  double Ee_th_brem;
  double EnergyIR;
  double Tempopt, Energyopt;
  /*Bremsstralhung */
  double Bremssxs;
  double totalsync, expectedsync;
  double totalIC,expectedIC;
  double ratio;
  if(source->continuous!=NULL){
    density=source->continuous->density;
    if (source->continuous->pulsartype==NULL)  Bfield=source->continuous->Bfield;
    else Bfield=getmagneticfield(source->continuous->time);
    TempIR=source->continuous->TempIR;
    EnergyIR=source->continuous->EnergyIR;
    Tempopt=source->continuous->Tempopt;
    Energyopt=source->continuous->Energyopt;
  }
  
  else if(source->impulsive!=NULL){
    density=source->impulsive->density;
    Bfield=source->impulsive->Bfield;
    TempIR=source->impulsive->TempIR;
    EnergyIR=source->impulsive->EnergyIR;
    Tempopt=source->impulsive->Tempopt;
    Energyopt=source->impulsive->Energyopt;
 }
  printf("density=%.2f\n",density);
  printf("Bfield=%.3e\n",Bfield);
  /*Maximum of the planck distribution derivation*/
 
  /*First the synchrotron emission */
  
  sizegamma=(20+8)*100/STEP+1;
  syncgamma=malloc(sizeof(struct ARRAYSIZE));
  ICgamma=malloc(sizeof(struct ARRAYSIZE));
  ICRgamma=malloc(sizeof(struct ARRAYSIZE));
  OPTgamma=malloc(sizeof(struct ARRAYSIZE));
  bremssgamma=malloc(sizeof(struct ARRAYSIZE));
  
  syncgamma->size=sizegamma;
  ICgamma->size=sizegamma;
  ICRgamma->size=sizegamma;
  OPTgamma->size=sizegamma;
  bremssgamma->size=sizegamma;
  for(i=0;i<2;i++){
    syncgamma->array[i]=malloc(sizegamma*sizeof(double));
    ICgamma->array[i]=malloc(sizegamma*sizeof(double));
    ICRgamma->array[i]=malloc(sizegamma*sizeof(double));
    OPTgamma->array[i]=malloc(sizegamma*sizeof(double));
    bremssgamma->array[i]=malloc(sizegamma*sizeof(double));
  }
  for(i=0;i<2;i++){
    for(j=0;j<sizegamma;j++){
      syncgamma->array[i][j]=0;
      ICgamma->array[i][j]=0;
      ICRgamma->array[i][j]=0;
      OPTgamma->array[i][j]=0;
      bremssgamma->array[i][j]=0;
    }
  }
  
for(i=0;i<=sizegamma;i++){
    syncgamma->array[0][i]=pow(10,(double)(i*STEP-2000.0)/100.0);
    ICgamma->array[0][i]=syncgamma->array[0][i];
    ICRgamma->array[0][i]=syncgamma->array[0][i];
    OPTgamma->array[0][i]=syncgamma->array[0][i];
    bremssgamma->array[0][i]=syncgamma->array[0][i];
    dEgamma=syncgamma->array[0][i]*(pow(10,(double)STEP/100.0)-1);
    //printf("ICgammaE=%.3e\n",ICgamma->array[0][i]);
    
    
    
    for(j=0;j<electronarraysize;j++){
      Ee=electrondistribution[0][j]*electronmassTeV;
      dEe=Ee*(pow(10,(double)STEP/100.0)-1); //separation between the 2 energy in array
      
      //printf("Ee=%.3e DEe=%.3e\n",Ee,dEe);
      syncgamma->array[1][i] +=Synchrotronphotonisotropic(Ee,syncgamma->array[0][i],Bfield)*electrondistribution[1][j]*dEe;

      /* Now for the IC distribution, we need the black body distribution*/
      /*We use the integration from Meyer et al to derive the one proton IC photon distribution*/
     
      
        if((Ee>ICgamma->array[0][i])&&(ICRgamma->array[0][i]>1E-12)){
	  IC1=0.0;
	  ICR=0.0;
	  ICO=0.0;
	Energyinitphoton=ICRgamma->array[0][i]/(4*pow(Ee/electronmassTeV,2))/(1-ICRgamma->array[0][i]/Ee);
	 

	 while (Energyinitphoton<ICRgamma->array[0][i]){
	dEnergyinit=Energyinitphoton*(pow(10,1.0/100.0)-1);
	//	printf("Energyinit=%.3e\n",Energyinitphoton);
	if (dEnergyinit>ICRgamma->array[0][i]-Energyinitphoton) dEnergyinit=ICgamma->array[0][i]-Energyinitphoton;
	/*for the CMB IC*/
       if(Energyinitphoton>ICgamma->array[0][i]/(4*pow(Ee/electronmassTeV,2))/(1-ICRgamma->array[0][i]/Ee)){
	  IC1 +=blackbodydistribution(Energyinitphoton,TEMPCMB)/Energyinitphoton*inversecomptonfunction(Ee, ICgamma->array[0][i], Energyinitphoton)*dEnergyinit;
	  //printf("IC1\n");
	  
	 // IC1 +=blackbodydistribution2(UCMBTeV*1E12,Energyinitphoton,TEMPCMB)/Energyinitphoton*inversecomptonfunction(Ee, ICgamma->array[0][i], Energyinitphoton)*dEnergyinit;
	}
	/*for the IR IC*/
	//printf("ICR\n");
	//printf("flush\n");
	if(Energyinitphoton>ICRgamma->array[0][i]/(4*pow(Ee/electronmassTeV,2))/(1-ICRgamma->array[0][i]/Ee)){
	ICR +=blackbodydistribution2(EnergyIR,Energyinitphoton,TempIR)/Energyinitphoton*inversecomptonfunction(Ee, ICgamma->array[0][i], Energyinitphoton)*dEnergyinit;
	}
	if(Energyinitphoton>OPTgamma->array[0][i]/(4*pow(Ee/electronmassTeV,2))/(1-OPTgamma->array[0][i]/Ee)){
	ICO +=blackbodydistribution2(Energyopt,Energyinitphoton,Tempopt)/Energyinitphoton*inversecomptonfunction(Ee, OPTgamma->array[0][i], Energyinitphoton)*dEnergyinit;
	}
	/*if(2.8*BOLTZMANN*TempIR/TeVtoJoules>ICgamma->array[0][i]/(4*pow(Ee/electronmassTeV,2))/(1-ICRgamma->array[0][i]/Ee)){
	  ICR=EnergyIR/pow(2.8*BOLTZMANN*TempIR/TeVtoJoules,2.0)*inversecomptonfunction(Ee,ICgamma->array[0][i],2.8*BOLTZMANN*TempIR/TeVtoJoules);
	  }*/

	
	Energyinitphoton +=dEnergyinit;
	 }
      
	 ICgamma->array[1][i] +=3.0/(4.0)*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*IC1*electrondistribution[1][j]*dEe;
	  ICRgamma->array[1][i] +=3.0/(4.0)*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*ICR*electrondistribution[1][j]*dEe;
	  OPTgamma->array[1][i] +=3.0/(4.0)*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*ICO*electrondistribution[1][j]*dEe;
	
	}
      Ee_th_brem=bremssgamma->array[0][i]+electronmassTeV;
      /*Now for the Bremsstralhung*/
      if (Ee > Ee_th_brem && bremssgamma->array[0][i] > 1.0e-15) {
	//Bremssxs = 0.75*BremssXS(bremssgamma->array[0][i], Ee, 1, 1) + 0.25*BremssXS(bremssgamma->array[0][i], Ee, 2, 2);
	Bremssxs =0.75*BremssXS(bremssgamma->array[0][i], Ee, 1, 1);
	bremssgamma->array[1][i] +=density*lightspeedcm*electrondistribution[1][j]*Bremssxs*dEe;
	
      }
      
    }
    totalsync +=1.0/ergtoTeV*syncgamma->array[1][i]*pow(syncgamma->array[0][i],1.0)*dEgamma;
    totalIC +=1.0/ergtoTeV*ICgamma->array[1][i]*pow(ICgamma->array[0][i],1.0)*dEgamma;

 }
 printf("totalIC= %.3e totalsync=%.3e ratio=%.3f expected=%.3e \n",totalIC,totalsync,totalIC/totalsync,UCMBTeV/ergtoTeV/pow(Bfield*1E-6,2.0)*8*PI);
  


}


void electronSEDsecondary(){
  /*The main difference between this function and the previous one is the fact the secondary electron energy is stored inside the array*/
 int i,j,l;
 double Egamma;
 double dEe,Ee;

 double density;
 double Bfield;
 double Energyinitphoton;
 double Ee_th_brem;
 double Bremssxs;
if(source->continuous!=NULL){
    density=source->continuous->density;
  if (source->continuous->pulsartype==NULL)  Bfield=source->continuous->Bfield;
  else Bfield=getmagneticfield(source->continuous->time);
 }

 else if(source->impulsive!=NULL){
    density=source->impulsive->density;
    Bfield=source->impulsive->Bfield;
 }


/*initialise secondaryelecgamma*/
 secondaryelecgammasync=malloc(sizeof(struct ARRAYSIZE));
 secondaryelecgammasync->size=syncgamma->size;
 secondaryelecgammabremss=malloc(sizeof(struct ARRAYSIZE));
 secondaryelecgammabremss->size=syncgamma->size;

for(i=0;i<2;i++){
    secondaryelecgammasync->array[i]=malloc(syncgamma->size*sizeof(double));
    secondaryelecgammabremss->array[i]=malloc(syncgamma->size*sizeof(double));
    for(j=0;j<syncgamma->size;j++){
      secondaryelecgammasync->array[i][j]=0.0;
      secondaryelecgammabremss->array[i][j]=0.0;
  }
 }
 
 for(i=0;i<syncgamma->size;i++){
  
 
  secondaryelecgammasync->array[0][i]=syncgamma->array[0][i];
  secondaryelecgammabremss->array[0][i]=syncgamma->array[0][i];
  printf("secondaryelecgamma=%.3e\n",secondaryelecgammasync->array[0][i]);
  /*Now let s use the secondary electron array secondaryelecgammasync*/
  for(j=0;j<secondaryelecpp->size;j++){
    
    Ee=secondaryelecpp->array[0][j];
    dEe=Ee*(pow(10,(double)STEP/100.0)-1);
    
    secondaryelecgammasync->array[1][i] +=Synchrotronphotonisotropic(Ee,secondaryelecgammasync->array[0][i],Bfield)*secondaryelecpp->array[1][j]*dEe;
    //printf("another useless quote secelecppsize=%d !\n",syncgamma->size);
      /* Now for the IC distribution, we need the black body distribution*/
      /*We use the integration from Meyer et al to derive the one proton IC photon distribution*/
      
      
 
      /*Now for the Bremsstralhung*/
    Ee_th_brem=secondaryelecgammasync->array[0][i]+electronmassTeV;
      if (Ee > Ee_th_brem && bremssgamma->array[0][i] > 1.0e-15) {
	Bremssxs = 0.75*BremssXS(bremssgamma->array[0][i], Ee, 1, 1);
     
	secondaryelecgammabremss->array[1][i] +=density*lightspeedcm*secondaryelecpp->array[1][i]*Bremssxs*dEe;

      }
  }
   
 }
  

}





double Synchrotronphotonisotropic(double Ee, double Egamma,double B){

  double criticalfrequency;
  double frequency;
 
  double Bessel,x;
  double fluxa;
  double totalflux=0.0;
  double alpha,dalpha=0.01;

  for(alpha=0.01;alpha<PI/2.0;alpha +=dalpha){
  frequency=Egamma*TeVtoJoules/PLANCK;
  criticalfrequency=3.0/2.0*echarge*B*microGausstoTesla*pow(lightspeedSI,2)/(2*PI*electronmassTeV*TeVtoJoules)*pow(Ee/electronmassTeV,2.0)*sin(alpha);
  
  
  //criticalfrequency=3.0/2.0*echarge*B*microGausstoTesla*pow(lightspeedSI,2)*sqrt(2.0/3.0)/(2*PI*electronmassTeV*TeVtoJoules)*pow(Ee/electronmassTeV,2.0);
  x=(double)frequency/criticalfrequency;
  //printf("sync limit= %.3e\n",-8.0*GSL_LOG_DBL_MIN/7.0);
  if ((x >= 0.0) && (x<-8.0*GSL_LOG_DBL_MIN/7.0)) { 
 
    Bessel=gsl_sf_synchrotron_1(x);
    
    if (Bessel<1E-9) Bessel=0.0;
  }
  else Bessel=0.0;
   


  //fluxa=sqrt(3.0)*pow(echarge,3.0)*B*sqrt(2.0/3.0)*microGausstoTesla*lightspeedSI/(16*pow(PI,2)*vacuumpermittivity*electronmassTeV*Egamma*pow(TeVtoJoules,2)*PLANCK); //ph s-1 J-1

  fluxa=sqrt(3.0)*pow(echarge,3.0)*B*microGausstoTesla*lightspeedSI/(4.0*pow(PI,1)*vacuumpermittivity*electronmassTeV*Egamma*pow(TeVtoJoules,2)*PLANCK)*sin(alpha); //ph s-1 J-1
  totalflux +=(double)(fluxa*Bessel)*TeVtoJoules*sin(alpha)*dalpha; //ph s-1 TeV-1
  }
  //  totalflux *=2.0; //integration over PI
  return totalflux;
}




double blackbodydistribution(double Egamma, double T){
  //returns number of the planckian distribution function describing the number of photons between Egamma and Egamma+dEgamma

  double distribution;
  //printf("Egamma=%.3e \n",Egamma);
  distribution=8.0*PI*pow(Egamma*TeVtoJoules,2)/(pow(PLANCK*lightspeedcm,3.0))*1.0/(exp(Egamma*TeVtoJoules/(BOLTZMANN*T))-1.0); //in ph cm-3 J-1
 
  //if (distribution<0)printf("distribution =%.3e\n",distribution);

  
  distribution *=TeVtoJoules; //in ph cm-3 TeV-1 
 
  return distribution;
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
  /*if (Egamma>1E-2 &&Egamma<1E-1) {
    printf("inverse compton=%.3f q=%.3f p =%.3f epsolon=%.3e\n",inversecompton,q,p,initphotonenergy);
    sleep(1);
    }*/
  if (q>1){

    // printf("q=%.3e Egamma=%.3e Ee=%.3e initphoton=%.3e\n",q,Egamma,Ee,initphotonenergy);
    inversecompton=0.0;
  }
  if(inversecompton<0) printf("IMPOSSIIIIIIIIIBLEEEEE\n");
  return inversecompton;

}


double ICcrosssectionthompson(double Ee,double Egamma,double energyinit){
  double inversecompton;
  inversecompton=2*Egamma*log(Egamma*pow(electronmassTeV/Ee,2)/(4*energyinit))+Egamma+4*pow(Ee/electronmassTeV,2)*energyinit-pow(Egamma,2)/(2*pow(Ee/electronmassTeV,2)*energyinit);


  return inversecompton;
}




/*ALL THAT REMAIN IS THE DIRAC FUNCTION FOR THE PRIMARY ELECTRON */
void electronSEDdiractest1(){
  /*This function has the same structure as electronSEDprimary with  the exception that thehre is no electrondistribution array as the it is monoenergetic*/
  int i,j,l;
  double Egamma,dEgamma;
  double Ee,dEe;
  
  double syncsize1,syncsize2;
  int sizegamma;
  
  double Bfield;
  double TempIR,EnergyIR;
  double Tempopt,Energyopt;
  double density;
  double dEnergyinit; //in TeV
  double EnergyinitphotonCMB,Energyinitphoton;
  double IC1,ICR,ICO;
  double Ee_th_brem;

  /*Bremsstralhung */
  double Bremssxs;

    density=source->dirac->density;
    Bfield=source->dirac->Bfield;
    TempIR=source->dirac->TempIR;
    EnergyIR=source->dirac->EnergyIR; //get the value in TeV
    Tempopt=source->dirac->Tempopt;
    Energyopt=source->dirac->Energyopt;
 

  /*Maximum of the planck distribution derivation*/
 
  /*First the synchrotron emission */
 
  sizegamma=(20+8)*100/STEP+1;
  syncgamma=malloc(sizeof(struct ARRAYSIZE));
  ICgamma=malloc(sizeof(struct ARRAYSIZE));
  ICRgamma=malloc(sizeof(struct ARRAYSIZE));
  OPTgamma=malloc(sizeof(struct ARRAYSIZE));
  bremssgamma=malloc(sizeof(struct ARRAYSIZE));
  
  syncgamma->size=sizegamma;
  ICgamma->size=sizegamma;
  ICRgamma->size=sizegamma;
  OPTgamma->size=sizegamma;
  bremssgamma->size=sizegamma;
  for(i=0;i<2;i++){
    syncgamma->array[i]=malloc(sizegamma*sizeof(double));
    ICgamma->array[i]=malloc(sizegamma*sizeof(double));
    ICRgamma->array[i]=malloc(sizegamma*sizeof(double));
    OPTgamma->array[i]=malloc(sizegamma*sizeof(double));
    bremssgamma->array[i]=malloc(sizegamma*sizeof(double));
  }
  for(i=0;i<2;i++){
    for(j=0;j<sizegamma;j++){
      syncgamma->array[i][j]=0;
      ICgamma->array[i][j]=0;
      ICRgamma->array[i][j]=0;
      OPTgamma->array[i][j]=0;
      bremssgamma->array[i][j]=0;
    }
  }
 
for(i=0;i<=sizegamma;i++){
    syncgamma->array[0][i]=pow(10,(double)(i*STEP-2000.0)/100.0);
    ICgamma->array[0][i]=syncgamma->array[0][i];
    ICRgamma->array[0][i]=syncgamma->array[0][i];
    OPTgamma->array[0][i]=syncgamma->array[0][i];
    bremssgamma->array[0][i]=syncgamma->array[0][i];
    dEgamma=gammaraypp->array[0][i]*(pow(10,(double)STEP/100.0)-1);
    
   
    
      Ee=electrondistribution[0][0]*electronmassTeV;
      // printf("Energy Ee=%.3e\n",Ee);
      dEe=Ee*(pow(10,(double)STEP/100.0)-1); //separation between the 2 energy in array
     
      syncgamma->array[1][i] +=Synchrotronphotonisotropic(Ee,syncgamma->array[0][i],Bfield)*electrondistribution[1][0];

      /* Now for the IC distribution, we need the black body distribution*/
      /*We use the integration from Meyer et al to derive the one proton IC photon distribution*/
      
      
      IC1=0.0;
      ICR=0.0;
      ICO=0.0;
      if((Ee>ICgamma->array[0][i])&&(ICRgamma->array[0][i]>1E-12)){

	Energyinitphoton=2.8*BOLTZMANN*2.73/TeVtoJoules/1000.0;
      /*for the CMB photon*/

	while (Energyinitphoton<2.8*BOLTZMANN*2.73/TeVtoJoules*10){
	 
	  dEnergyinit=2.8*BOLTZMANN*2.73/TeVtoJoules/100.0;
	 
	/*for the CMB IC*/
	   if(Energyinitphoton>ICRgamma->array[0][i]/(4*pow(Ee/electronmassTeV,2))/(1-ICRgamma->array[0][i]/Ee)){
	    
	     IC1 +=blackbodydistribution2(UCMBTeV*1E12,Energyinitphoton,TEMPCMB)/(Energyinitphoton)*inversecomptonfunction(Ee, syncgamma->array[0][i], Energyinitphoton)*dEnergyinit;
	     //IC1 +=blackbodydistribution2(UCMBTeV*1E12,Energyinitphoton,TEMPCMB)/pow(Energyinitphoton,2)*ICcrosssectionthompson(Ee, syncgamma->array[0][i], Energyinitphoton)*dEnergyinit;
	}
	   Energyinitphoton +=dEnergyinit;
      }
	//printf("IC=%.3e\n",IC1);
	
	/*for the IR IC*/
	Energyinitphoton=2.8*BOLTZMANN*TempIR/TeVtoJoules/1000.0;
      
	//	if((ICgamma->array[0][i]>20)&&(ICgamma->array[0][i]<30)) sleep(5);
	while (Energyinitphoton<2.8*BOLTZMANN*TempIR/TeVtoJoules*10){
	 
	  dEnergyinit=2.8*BOLTZMANN*TempIR/TeVtoJoules/100.0;
	  
	/*for the CMB IC*/

	  if(Energyinitphoton>ICRgamma->array[0][i]/(4*pow(Ee/electronmassTeV,2))/(1-ICRgamma->array[0][i]/Ee)){
	    
	    ICR +=blackbodydistribution2(EnergyIR,Energyinitphoton,TempIR)/Energyinitphoton*inversecomptonfunction(Ee, ICgamma->array[0][i], Energyinitphoton)*dEnergyinit;
	  }	
	  Energyinitphoton +=dEnergyinit;
	}

	/*for the optical starlight */
		while (Energyinitphoton<2.8*BOLTZMANN*Tempopt/TeVtoJoules*10){
	 
	  dEnergyinit=2.8*BOLTZMANN*TempIR/TeVtoJoules/100.0;
	  
	/*for the CMB IC*/

	  if(Energyinitphoton>OPTgamma->array[0][i]/(4*pow(Ee/electronmassTeV,2))/(1-OPTgamma->array[0][i]/Ee)){
	    
	    ICO +=blackbodydistribution2(Energyopt,Energyinitphoton,Tempopt)/Energyinitphoton*inversecomptonfunction(Ee, OPTgamma->array[0][i], Energyinitphoton)*dEnergyinit;
	  }	
	  Energyinitphoton +=dEnergyinit;
	}

	ICgamma->array[1][i]=3.0/4.0*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*IC1*electrondistribution[1][0];
	//ICgamma->array[1][i]=3.0/16.0*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,4)*IC1*electrondistribution[1][0];
	ICRgamma->array[1][i]=3.0/4.0*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*ICR*electrondistribution[1][0];
	OPTgamma->array[1][i]=3.0/4.0*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*ICO*electrondistribution[1][0];
      
      }
      //if((ICgamma->array[0][i]>20)&&(ICgamma->array[0][i]<30)) sleep(5);
      Ee_th_brem=bremssgamma->array[0][i]+electronmassTeV;
      /*Now for the Bremsstralhung*/
      if (Ee > Ee_th_brem && bremssgamma->array[0][i] > 1.0e-15) {
	Bremssxs =BremssXS(bremssgamma->array[0][i], Ee, 1, 1);
     
	bremssgamma->array[1][i] +=density*lightspeedcm*electrondistribution[1][0]*Bremssxs;

      }
      
    

}

  


}






void electronSEDdirac(){
  /*This function has the same structure as electronSEDprimary with  the exception that thehre is no electrondistribution array as the it is monoenergetic*/
  int i,j,l;
  double Egamma,dEgamma;
  double Ee,dEe;
  
  double syncsize1,syncsize2;
  int sizegamma;
  
  double Bfield;
  double TempIR,EnergyIR;
   double Tempopt,Energyopt;
  double density;
  double dEnergyinit; //in TeV
  double EnergyinitphotonCMB,Energyinitphoton;
  double IC1,ICR,ICO;
  double Ee_th_brem;

  /*Bremsstralhung */
  double Bremssxs;
  double totalIC,totalsync;
  double expectedsync,expectedIC;

    density=source->dirac->density;
    Bfield=source->dirac->Bfield;
    TempIR=source->dirac->TempIR;
    EnergyIR=source->dirac->EnergyIR; //get the value in TeV
    Tempopt=source->dirac->Tempopt;
    Energyopt=source->dirac->Energyopt;
 

  /*Maximum of the planck distribution derivation*/
 
  /*First the synchrotron emission */
 
  sizegamma=(20+8)*100/STEP+1;
  syncgamma=malloc(sizeof(struct ARRAYSIZE));
  ICgamma=malloc(sizeof(struct ARRAYSIZE));
  ICRgamma=malloc(sizeof(struct ARRAYSIZE));
  OPTgamma=malloc(sizeof(struct ARRAYSIZE));
  bremssgamma=malloc(sizeof(struct ARRAYSIZE));
  
  syncgamma->size=sizegamma;
  ICgamma->size=sizegamma;
  ICRgamma->size=sizegamma;
  OPTgamma->size=sizegamma;
  bremssgamma->size=sizegamma;
  for(i=0;i<2;i++){
    syncgamma->array[i]=malloc(sizegamma*sizeof(double));
    ICgamma->array[i]=malloc(sizegamma*sizeof(double));
    ICRgamma->array[i]=malloc(sizegamma*sizeof(double));
    OPTgamma->array[i]=malloc(sizegamma*sizeof(double));
    bremssgamma->array[i]=malloc(sizegamma*sizeof(double));
  }
  for(i=0;i<2;i++){
    for(j=0;j<sizegamma;j++){
      syncgamma->array[i][j]=0;
      ICgamma->array[i][j]=0;
      ICRgamma->array[i][j]=0;
      OPTgamma->array[i][j]=0;
      bremssgamma->array[i][j]=0;
    }
  }
 
for(i=0;i<=sizegamma;i++){
    syncgamma->array[0][i]=pow(10,(double)(i*STEP-2000.0)/100.0);
    ICgamma->array[0][i]=syncgamma->array[0][i];
    ICRgamma->array[0][i]=syncgamma->array[0][i];
    OPTgamma->array[0][i]=syncgamma->array[0][i];
    bremssgamma->array[0][i]=syncgamma->array[0][i];
    dEgamma=gammaraypp->array[0][i]*(pow(10,(double)STEP/100.0)-1);
    
   
    
      Ee=electrondistribution[0][0]*electronmassTeV;
     
      // printf("Energy Ee=%.3e\n",Ee);
      dEe=Ee*(pow(10,(double)STEP/100.0)-1); //separation between the 2 energy in array
     
      syncgamma->array[1][i] +=Synchrotronphotonisotropic(Ee,syncgamma->array[0][i],Bfield)*electrondistribution[1][0];

      /* Now for the IC distribution, we need the black body distribution*/
      /*We use the integration from Meyer et al to derive the one proton IC photon distribution*/
      totalsync +=1.0/ergtoTeV*syncgamma->array[1][i]*pow(syncgamma->array[0][i],1.0)*dEgamma;
      
      IC1=0.0;
      ICR=0.0;
      ICO=0.0;
      if((Ee>ICgamma->array[0][i])&&(ICRgamma->array[0][i]>1E-12)){
	 IC1=0.0;
	  ICR=0.0;
	  ICO=0.0;
	Energyinitphoton=ICRgamma->array[0][i]/(4*pow(Ee/electronmassTeV,2))/(1-ICRgamma->array[0][i]/Ee);
	 

	 while (Energyinitphoton<ICRgamma->array[0][i]){
	dEnergyinit=Energyinitphoton*(pow(10,1.0/100.0)-1);
	//	printf("Energyinit=%.3e\n",Energyinitphoton);
	if (dEnergyinit>ICRgamma->array[0][i]-Energyinitphoton) dEnergyinit=ICgamma->array[0][i]-Energyinitphoton;
	/*for the CMB IC*/
       if(Energyinitphoton>ICgamma->array[0][i]/(4*pow(Ee/electronmassTeV,2))/(1-ICRgamma->array[0][i]/Ee)){
	  
	  //printf("IC1\n");
	  IC1 +=blackbodydistribution2(UCMBTeV*1E12,Energyinitphoton,TEMPCMB)/pow(Energyinitphoton,2)*ICcrosssectionthompson(Ee, syncgamma->array[0][i], Energyinitphoton)*dEnergyinit;
	 // IC1 +=blackbodydistribution2(UCMBTeV*1E12,Energyinitphoton,TEMPCMB)/Energyinitphoton*inversecomptonfunction(Ee, ICgamma->array[0][i], Energyinitphoton)*dEnergyinit;
	}
	/*for the IR IC*/
	//printf("ICR\n");
	//printf("flush\n");
	if(Energyinitphoton>ICRgamma->array[0][i]/(4*pow(Ee/electronmassTeV,2))/(1-ICRgamma->array[0][i]/Ee)){
	ICR +=blackbodydistribution2(EnergyIR,Energyinitphoton,TempIR)/Energyinitphoton*inversecomptonfunction(Ee, ICgamma->array[0][i], Energyinitphoton)*dEnergyinit;
	}
	if(Energyinitphoton>OPTgamma->array[0][i]/(4*pow(Ee/electronmassTeV,2))/(1-OPTgamma->array[0][i]/Ee)){
	ICO +=blackbodydistribution2(EnergyIR,Energyinitphoton,TempIR)/Energyinitphoton*inversecomptonfunction(Ee, ICgamma->array[0][i], Energyinitphoton)*dEnergyinit;
	}
	/*if(2.8*BOLTZMANN*TempIR/TeVtoJoules>ICgamma->array[0][i]/(4*pow(Ee/electronmassTeV,2))/(1-ICRgamma->array[0][i]/Ee)){
	  ICR=EnergyIR/pow(2.8*BOLTZMANN*TempIR/TeVtoJoules,2.0)*inversecomptonfunction(Ee,ICgamma->array[0][i],2.8*BOLTZMANN*TempIR/TeVtoJoules);
	  }*/

	
	Energyinitphoton +=dEnergyinit;
	 }
	
	  ICgamma->array[1][i] =3.0/(4.0)*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*IC1*electrondistribution[1][0];
	 //ICgamma->array[1][i]=3.0/16.0/(4*PI)*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,4)*IC1*electrondistribution[1][0];
	 ICRgamma->array[1][i] =3.0/(4.0)*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*ICR*electrondistribution[1][0];	
	 OPTgamma->array[1][i] =3.0/(4.0)*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*ICO*electrondistribution[1][0];	
      }
      totalIC +=1.0/ergtoTeV*ICgamma->array[1][i]*pow(ICgamma->array[0][i],1.0)*dEgamma;
	//printf("IC=%.3e\n",IC1);
	
	/*for the IR IC*/

      //if((ICgamma->array[0][i]>20)&&(ICgamma->array[0][i]<30)) sleep(5);
      Ee_th_brem=bremssgamma->array[0][i]+electronmassTeV;
      /*Now for the Bremsstralhung*/
      if (Ee > Ee_th_brem && bremssgamma->array[0][i] > 1.0e-15) {
	Bremssxs =0.75*BremssXS(bremssgamma->array[0][i], Ee, 1, 1);
     
	bremssgamma->array[1][i] +=4*PI*density*lightspeedcm*electrondistribution[1][0]*Bremssxs;

      }
      
    

}

 
 expectedsync=4.0/3.0*thompsonXS*lightspeedcm*pow(source->dirac->Edirac/electronmassTeV,2.0)*pow(Bfield*1E-6,2.0)/(8*PI)*electrondistribution[1][0];
 expectedIC=1.0/ergtoTeV*4.0/3.0*thompsonXS*lightspeedcm*pow(source->dirac->Edirac/electronmassTeV,2.0)*UCMBTeV*electrondistribution[1][0];;
 printf(" totalICenergy=%.3e ------totalenergysync=%.3e \n ratio=%.3e expected ratio= %.3e \n",totalIC,totalsync,totalIC/totalsync,UCMBTeV/ergtoTeV*8*PI/pow(Bfield*1E-6,2));
 printf("expectedIC=%.3e -----expectedsync=%.3e ratio=%.3f\n",expectedIC, expectedsync,expectedIC/expectedsync);
 printf("totalerrorIC=%.3f totalerrorsync=%.3f\n",(4*PI)*totalIC/expectedIC,(4*PI)*totalsync/expectedsync);
}







void electronSEDdiractest2(){
  /*This function has the same structure as electronSEDprimary with  the exception that thehre is no electrondistribution array as the it is monoenergetic*/
  int i,j,l;
  double Egamma,dEgamma;
  double Ee,dEe;
  
  double syncsize1,syncsize2;
  int sizegamma;
  
  double Bfield;
  double TempIR,EnergyIR;
  double density;
  double dEnergyinit; //in TeV
  double EnergyinitphotonCMB,Energyinitphoton;
  double IC1,ICR;
  double Ee_th_brem;

  /*Bremsstralhung */
  double Bremssxs;
  double totalIC=0.0;
  double totalsync=0.0;
    density=source->dirac->density;
    Bfield=source->dirac->Bfield;
    TempIR=source->dirac->TempIR;
    EnergyIR=source->dirac->EnergyIR; //get the value in TeV
 
 

  /*Maximum of the planck distribution derivation*/
 
  /*First the synchrotron emission */
 
  sizegamma=(20+8)*100/STEP+1;
  syncgamma=malloc(sizeof(struct ARRAYSIZE));
  ICgamma=malloc(sizeof(struct ARRAYSIZE));
  ICRgamma=malloc(sizeof(struct ARRAYSIZE));
  bremssgamma=malloc(sizeof(struct ARRAYSIZE));
  
  syncgamma->size=sizegamma;
  ICgamma->size=sizegamma;
  ICRgamma->size=sizegamma;
  bremssgamma->size=sizegamma;
  for(i=0;i<2;i++){
    syncgamma->array[i]=malloc(sizegamma*sizeof(double));
    ICgamma->array[i]=malloc(sizegamma*sizeof(double));
    ICRgamma->array[i]=malloc(sizegamma*sizeof(double));
    bremssgamma->array[i]=malloc(sizegamma*sizeof(double));
  }
  for(i=0;i<2;i++){
    for(j=0;j<sizegamma;j++){
      syncgamma->array[i][j]=0;
      ICgamma->array[i][j]=0;
      ICRgamma->array[i][j]=0;
      bremssgamma->array[i][j]=0;
    }
  }
 
for(i=0;i<=sizegamma;i++){
    syncgamma->array[0][i]=pow(10,(double)(i*STEP-2000.0)/100.0);
    ICgamma->array[0][i]=syncgamma->array[0][i];
    ICRgamma->array[0][i]=syncgamma->array[0][i];
    bremssgamma->array[0][i]=syncgamma->array[0][i];
    dEgamma=gammaraypp->array[0][i]*(pow(10,(double)STEP/100.0)-1);
    
   
    
      Ee=electrondistribution[0][0]*electronmassTeV;
      // printf("Energy Ee=%.3e\n",Ee);
      dEe=Ee*(pow(10,(double)STEP/100.0)-1); //separation between the 2 energy in array
     
      syncgamma->array[1][i] +=Synchrotronphotonisotropic(Ee,syncgamma->array[0][i],Bfield)*electrondistribution[1][0];
      totalsync +=1.0/ergtoTeV*syncgamma->array[1][i]*pow(syncgamma->array[0][i],2.0)*dEgamma;
      /* Now for the IC distribution, we need the black body distribution*/
      /*We use the integration from Meyer et al to derive the one proton IC photon distribution*/
      
      
      
      IC1=0.0;
      ICR=0.0;
      Energyinitphoton=ICRgamma->array[0][i]/(4*pow(Ee/electronmassTeV,2))/(1-ICRgamma->array[0][i]/Ee);
	 

	
	//	printf("Energyinit=%.3e\n",Energyinitphoton);

	/*for the CMB IC*/
    
	
	  
      IC1 =UCMBTeV/pow(2.8*BOLTZMANN*TEMPCMB/TeVtoJoules,2)*inversecomptonfunction(Ee, ICgamma->array[0][i], 2.8*BOLTZMANN*TEMPCMB/TeVtoJoules);
     
	/*for the IR IC*/
      ICR =EnergyIR/pow(2.8*BOLTZMANN*TempIR/TeVtoJoules,2)*inversecomptonfunction(Ee, ICgamma->array[0][i], 2.8*BOLTZMANN*TempIR/TeVtoJoules);

 
	/*if(2.8*BOLTZMANN*TempIR/TeVtoJoules>ICgamma->array[0][i]/(4*pow(Ee/electronmassTeV,2))/(1-ICRgamma->array[0][i]/Ee)){
	  ICR=EnergyIR/pow(2.8*BOLTZMANN*TempIR/TeVtoJoules,2.0)*inversecomptonfunction(Ee,ICgamma->array[0][i],2.8*BOLTZMANN*TempIR/TeVtoJoules);
	  }*/

	
      
 
      
	  ICgamma->array[1][i] =3.0/4.0*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*IC1*electrondistribution[1][0];
	 ICRgamma->array[1][i] =3.0/4.0*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*ICR*electrondistribution[1][0];	
	
	  totalIC +=1.0/ergtoTeV*ICgamma->array[1][i]*pow(ICgamma->array[0][i],2.0)*dEgamma;
	 ;
	//printf("IC=%.3e\n",IC1);
	
	/*for the IR IC*/

      //if((ICgamma->array[0][i]>20)&&(ICgamma->array[0][i]<30)) sleep(5);
      Ee_th_brem=bremssgamma->array[0][i]+electronmassTeV;
      /*Now for the Bremsstralhung*/
      if (Ee > Ee_th_brem && bremssgamma->array[0][i] > 1.0e-15) {
	Bremssxs =BremssXS(bremssgamma->array[0][i], Ee, 1, 1);
     
	bremssgamma->array[1][i] +=density*lightspeedcm*electrondistribution[1][0]*Bremssxs;

      }
      
    

 }

  


}




void SSCgammafunction(){
  double size;
  int i,j,l;
  double dEgamma;
  double tempSSC;
  double initphotonSSC;
  double pSSC;
  double dEnergyinit;
  double Ee,dEe;
  if(source->continuous!=NULL){
    size=source->continuous->size*parsectocm; //size in cm for relevance to density
  }

  else if(source->impulsive!=NULL){
    size=source->impulsive->size*parsectocm;
  }
  SSCgamma=malloc(sizeof(struct ARRAYSIZE));
  SSCgamma->size=syncgamma->size;
  for(i=0;i<2;i++){
    SSCgamma->array[i]=malloc(SSCgamma->size*sizeof(double));
    for(j=0;j<SSCgamma->size;j++){
      SSCgamma->array[i][j]=0.0;
    }
  }
  

  for (i=0;i<SSCgamma->size;i++){
   
    SSCgamma->array[0][i]=pow(10,(double)(i*STEP-2000.0)/100.0);
    dEgamma=SSCgamma->array[0][i]*(pow(10,(double)STEP/100.0)-1);
    
    for (j=0;j<SSCgamma->size;j++){
      	initphotonSSC=syncgamma->array[0][j];
	dEnergyinit=initphotonSSC*(pow(10,(double)STEP/100.0)-1);
      for(l=0;l<electronarraysize;l++){	
	Ee=electrondistribution[0][l]*electronmassTeV;
	dEe=Ee*(pow(10,(double)(STEP/100.0))-1);
      

	/*first the synchrotron primary*/
	tempSSC=3.0/4.0*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*inversecomptonfunction(Ee, SSCgamma->array[0][i], initphotonSSC)*electrondistribution[1][l]/initphotonSSC*syncgamma->array[1][j]/(4.0/3.0*PI*pow(size,3))*dEe;
	
	pSSC +=tempSSC;
	syncgamma->array[1][j]-=tempSSC;
	/*for the CMB IC*/
	tempSSC=3.0/4.0*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*inversecomptonfunction(Ee, SSCgamma->array[0][i], initphotonSSC)*electrondistribution[1][l]/initphotonSSC*ICgamma->array[1][j]/(4.0/3.0*PI*pow(size,3))*dEe;
	
	pSSC +=tempSSC;
	ICgamma->array[1][j] -=tempSSC;
    
	/*for the IR IC*/
	tempSSC=3.0/4.0*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*inversecomptonfunction(Ee, SSCgamma->array[0][i], initphotonSSC)*electrondistribution[1][l]/initphotonSSC*ICRgamma->array[1][j]/(4.0/3.0*PI*pow(size,3))*dEe;
	
	pSSC +=tempSSC;
	ICRgamma->array[1][j] -=tempSSC;
	
	  /*for the bremss*/
	  tempSSC=3.0/4.0*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*inversecomptonfunction(Ee, SSCgamma->array[0][i], initphotonSSC)*electrondistribution[1][l]/initphotonSSC*bremssgamma->array[1][j]/(4.0/3.0*PI*pow(size,3))*dEe;
	
	pSSC +=tempSSC;
	bremssgamma->array[1][j] -=tempSSC;

	  /*for secondary electron sync*/
	  
	 tempSSC=3.0/4.0*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*inversecomptonfunction(Ee, SSCgamma->array[0][i], initphotonSSC)*electrondistribution[1][l]/initphotonSSC*secondaryelecgammasync->array[1][j]/(4.0/3.0*PI*pow(size,3))*dEe;
	
	pSSC +=tempSSC;
	secondaryelecgammasync->array[1][j] -=tempSSC;

	/*for secondary electron sync*/
	  
	 tempSSC=3.0/4.0*thompsonXS*lightspeedcm*pow(electronmassTeV/Ee,2)*inversecomptonfunction(Ee, SSCgamma->array[0][i], initphotonSSC)*electrondistribution[1][l]/initphotonSSC*secondaryelecgammabremss->array[1][j]/(4.0/3.0*PI*pow(size,3))*dEe;
	
	pSSC +=tempSSC;
	secondaryelecgammabremss->array[1][j] -=tempSSC;

      }
      SSCgamma->array[1][i] +=pSSC*dEnergyinit;
    }
  }

}
