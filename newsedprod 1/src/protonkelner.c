#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structure.h"
#include "constant.h"
#include <string.h>
#include "globalparameters.h"
#include "protondistributionlib.h"
#include "protonSEDlib.h"

void ppinteractiongammakelner(){
  /*This functions aims at getting the main secondary electrons+neutrinos+gamma-ray from the given proton energy distribution, we remember that the dynamic array start with GMIN or gamma2 (see Manolakou et al 2007 for gamma2), thus we must expect */

  /*We are reminded that the main  energy unit here is TeV*/

  /*First the gammaloop*/
  int i,j; //loop integer variable
  double dEgamma;
  double dEp,Ep;
  double Tp,dTp;
  double x,dx;
  double *GMINMAX;
  double EGMIN;
  int sizegamma;
  /* for the secondary electron*/
  
  int sizesecelectron;
  
  double dEsec;

  /*for the neutrino*/
  int sizeneutrino;
 double dEneut;

  double density;
  if(source->continuous!=NULL) density=source->continuous->density;
  else if(source->impulsive!=NULL) density=source->impulsive->density;
  GMINMAX=protonGMINGMAX();
  EGMIN=GMINMAX[0]*protonmassTeV;
 
  /*initialisation of the gammaraypp*/
  sizegamma=(int)(20+8)*100/STEP;
  gammaraypp=malloc(sizeof(struct ARRAYSIZE));
  gammaraypp->size=sizegamma;
  for(i=0;i<2;i++){
    gammaraypp->array[i]=malloc(sizegamma*sizeof(double));
  }
  for(i=0;i<2;i++){
    for(j=0;j<sizegamma;j++){
      gammaraypp->array[i][j]=0;
    }
  }
  /*First the gamma-ray from proton-proton innteraction*/

  
  for(i=0;i<=sizegamma;i++){
    gammaraypp->array[0][i]=pow(10,(double)(i*STEP-2000.0)/100.0);

    dEgamma=gammaraypp->array[0][i]*(pow(10,(double)STEP/100.0)-1);
    
    
    
    for(j=0;j<protonarraysize;j++){
 
      Ep=protondistribution[0][j]*protonmassTeV;
      dEp=Ep*(pow(10,(double)STEP/100.0)-1); //separation between the 2 energy in array
      Tp=Ep-protonmassTeV;
      dTp=Tp*(pow(10,(double)STEP/100.0)-1);
      
      if ((Tp>0.2797E-3)&&(gammaraypp->array[0][i]<Ep)) gammaraypp->array[1][i] +=lightspeedcm*density*protondistribution[1][j]*PPCrossSection(Ep)*Fgamma(gammaraypp->array[0][i]/Ep,Ep)*dEp/Ep;
      

       

    }

  }

  
  sizesecelectron=(int)(20+8)*100/STEP+1;

  secondaryelecpp=malloc(sizeof(struct ARRAYSIZE));
  secondaryelecpp->size=sizesecelectron;
  for(i=0;i<2;i++){
    secondaryelecpp->array[i]=malloc(sizesecelectron*sizeof(double));
  }
  for(i=0;i<2;i++){
    for(j=0;j<sizegamma;j++){
      secondaryelecpp->array[i][j]=0;
    }
  }

  /*The secondary electron coming from proton-proton interaction*/

  
  for(i=0;i<=sizesecelectron;i++){
    secondaryelecpp->array[0][i]=pow(10,(double)(i*STEP-2000)/100.0);

    dEsec=secondaryelecpp->array[0][i]*(pow(10,(double)STEP/100.0)-1);
    
    
    
    
    for(j=0;j<protonarraysize;j++){
      Ep=protondistribution[0][j]*protonmassTeV;
      dEp=Ep*(pow(10,(double)STEP/100.0)-1); //separation between the 2 energy in array

   
      if((Ep>=secondaryelecpp->array[0][i])&&(secondaryelecpp->array[0][i]>electronmassTeV)){
	secondaryelecpp->array[1][i] +=lightspeedcm*density*protondistribution[1][j]*PPCrossSection(Ep)*Felectron(secondaryelecpp->array[0][i]/Ep,Ep)*dEp/Ep;

      }

    }

  }

  /*This part with the neutrino flux will need to be tested thoroughly, I need to look more thoroughly at the method needed to compute the different functions */
  sizeneutrino=28*100.0/STEP+1;
  neutrinospp=malloc(sizeof(struct ARRAYSIZE));
  neutrinospp->size=sizeneutrino;
  for(i=0;i<2;i++){
    neutrinospp->array[i]=malloc(sizeneutrino*sizeof(double));
  }
  for(i=0;i<2;i++){
    for(j=0;j<sizegamma;j++){
      neutrinospp->array[i][j]=0;
    }
  }





for(i=0;i<=sizeneutrino;i++){
    neutrinospp->array[0][i]=pow(10,(double)(i*STEP-8000.0)/100.0);

    dEneut=gammaraypp->array[0][i]*(pow(10,(double)STEP/100.0)-1);
    
    
    
    
    for(j=0;j<protonarraysize;j++){
      Ep=protondistribution[0][j]*protonmassTeV;
      dEp=Ep*(pow(10,(double)STEP/100.0)-1); //separation between the 2 energy in array
      if((Ep>=neutrinospp->array[0][i])&&(neutrinospp->array[0][i]>neutrinomassTeV)){
	neutrinospp->array[1][i] +=lightspeedcm*density*PPCrossSection(Ep)*protondistribution[1][j]*(Fgamma(neutrinospp->array[0][i]/Ep,Ep)+Fneutrino(neutrinospp->array[0][i]/Ep,Ep))*dEp/Ep;

      }

    }

 }

}


