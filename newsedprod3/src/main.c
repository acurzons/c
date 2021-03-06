#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
/*variables*/
#include "constant.h"
#include "structure.h"
#include "globalparameters.h"
/*functions*/
#include "coolingelectronlib.h"
#include "coolingprotonlib.h"
#include "electrondistributionlib.h"
#include "protondistributionlib.h"
#include "electronSEDlib.h"
#include "protonSEDlib.h"
#include "protonSEDlibdermer.h"
#include "getvaluelib.h"
#include "displayfunctionslib.h"
#include "pulsarmisclib.h"
int main(int argc , char **argv){
  /*We start by picking the different parameters and store them into our struct variable */
  char file[50];
  sprintf(file,"%s",argv[1]);
  /*there will be one inline argument for the moment which will be the prefix of the different file created in OUTPUT ---SEE DISPLAYfunctions---- */
  source=malloc(sizeof(struct sourcetype));
  getlistparameters();
  
   
  
  
printf("NEWSTEOre\n");
 
printf("NOw let s get to business\n");
 printf("now to electron\n");
 
  /*Case where the source is impulsive*/
  if(source->impulsive!=NULL){
    if(strcmp(source->impulsive->evolution,"no")==0){
      constantprotondistribution();
      constantelectrondistribution(); 
    }

    else if(strcmp(source->impulsive->evolution,"yes")==0){
      
      evolutionprotonimpulsive();
      evolutionelectronimpulsive();
    }
    /*Now get to the SED*/
    ppinteractiongamma();
    electronSEDprimary();
    electronSEDsecondary();
    if(SSCoption==1) SSCgammafunction();

  }
  
  else if(source->continuous!=NULL){
   if(source->continuous->pulsartype!=NULL){
    getBfieldlookuptable();
    } 
    if(strcmp(source->continuous->evolution,"no")==0){
      
      constantprotondistribution();
      constantelectrondistribution(); 
    }

    else if(strcmp(source->continuous->evolution,"yes")==0){
      
      evolutionprotoncontinuous();
     
      evolutionelectroncontinuous();
    }
    ppinteractiongamma();
    //ppinteractiongamma();
    electronSEDprimary();
    electronSEDsecondary();
    printf("get to the OUTPUT\n");
    if(SSCoption==1) SSCgammafunction();

 }

  else if(source->dirac!=NULL){
    
    diracprotondistribution();
    
    diracelectrondistribution();
    
    printf("We finished with the dirac distribution\n");
    ppinteractiondirac();
    electronSEDdirac();
  }

  makeoutput(file);
  /*Now onto the display*/
  return 0;
}
