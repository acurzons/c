#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "structure.h"
#include "globalparameters.h"
#include "getvaluelib.h"
#include "constant.h"
  
void getlistparameters(){
  FILE *tmp; //temporary file for 
  FILE *infile;
  char command[200];
  char filename[200];
  char protonspectrum[200];
  char electronspectrum[200];
  char pulsarname[200];
  
  sprintf(filename,"%sParameters/parameterlist.txt",getenv("SEDPRODFOLDER"));
  sprintf(electronspectrum,"%sParameters/electronspectra.txt",getenv("SEDPRODFOLDER"));
  sprintf(protonspectrum,"%sParameters/protonspectra.txt",getenv("SEDPRODFOLDER"));
  

  sprintf(pulsarname,"%sParameters/pulsarparam.txt",getenv("SEDPRODFOLDER"));
  /*this command will determine the type of source we are dealing with*/
  sprintf(command,"grep 'type' %s  | awk -F '=' '{print $2}' >tmp.txt",filename);
  system(command); 
if ((infile=fopen("tmp.txt","r"))==NULL){
    printf("Unable to read tmp.txt\n");
    exit(1);
  }

  fscanf(infile,"%s",source->sourcetype);
  printf("ola\n");
  printf("type %s\n",source->sourcetype);
  /*1st case : impulsive option*/
  if (strcmp(source->sourcetype,"impulsive")==0){
    /*allocate memory for the pointer*/
    source->impulsive=malloc(sizeof(struct impulsive));
    sprintf(command,"grep 'Energy' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'evolution' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'escape' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'Efractionproton' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'distance' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'size' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'density' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'Bfield' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'time' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'TempIR' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'EnIR' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
sprintf(command,"grep 'Tempopt' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'Enopt' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'spectratypeelectron' %s | awk -F '=' '{print $2}' >>tmp.txt",electronspectrum);
    system(command); 
    sprintf(command,"grep 'spectratypeproton' %s | awk -F '=' '{print $2}' >>tmp.txt",protonspectrum);
    system(command); 

    /*Now let's scan the different values*/
    fscanf(infile,"%lf",&(source->impulsive->E0));
    fscanf(infile,"%s",source->impulsive->evolution);
    fscanf(infile,"%s",source->impulsive->escape);
    fscanf(infile,"%lf",&(source->impulsive->Efractionproton));
    fscanf(infile,"%lf",&(source->impulsive->distance));
    fscanf(infile,"%lf",&(source->impulsive->size));
    fscanf(infile,"%lf",&(source->impulsive->density));
    fscanf(infile,"%lf",&(source->impulsive->Bfield));
    fscanf(infile,"%lf",&(source->impulsive->time));
    source->impulsive->time *=365.24*3600*24;
    fscanf(infile,"%lf",&(source->impulsive->TempIR));
    fscanf(infile,"%lf",&(source->impulsive->EnergyIR));
    fscanf(infile,"%lf",&(source->impulsive->Tempopt));
    fscanf(infile,"%lf",&(source->impulsive->Energyopt));
    fscanf(infile,"%s",source->impulsive->spectratypeelectron);
    fscanf(infile,"%s",source->impulsive->spectratypeproton);
    

    /*  electron spectra inut type */
    if (strcmp(source->impulsive->spectratypeelectron,"powerlaw")==0){
      source->impulsive->powerlawelectron=malloc(sizeof(struct powerlaw));
      source->impulsive->powerlawelectron=getspectrumpowerlaw(infile,electronspectrum);
      source->impulsive->powerlawelectron->GMIN=source->impulsive->powerlawelectron->EMIN/electronmassTeV +1.0;
      source->impulsive->powerlawelectron->GMAX=source->impulsive->powerlawelectron->EMAX/electronmassTeV +1.0;
      
    }
    
    else if(strcmp(source->impulsive->spectratypeelectron,"brokenpowerlaw")==0){
      source->impulsive->brokenpowerlawelectron=malloc(sizeof(struct brokenpowerlaw));
      source->impulsive->brokenpowerlawelectron=getspectrumbrokenpowerlaw(infile,electronspectrum);
      source->impulsive->brokenpowerlawelectron->GMIN=source->impulsive->brokenpowerlawelectron->EMIN/electronmassTeV +1.0;
      source->impulsive->brokenpowerlawelectron->GMAX=source->impulsive->brokenpowerlawelectron->EMAX/electronmassTeV +1.0;

    }
    
    else if(strcmp(source->impulsive->spectratypeelectron,"powerlawcutoff")==0){
     source->impulsive->powerlawcutoffelectron=malloc(sizeof(struct powerlawcutoff)); 
     source->impulsive->powerlawcutoffelectron=getspectrumpowerlawcutoff(infile,electronspectrum);
     source->impulsive->powerlawcutoffelectron->GMIN=source->impulsive->powerlawcutoffelectron->EMIN/electronmassTeV +1.0;
     source->impulsive->powerlawcutoffelectron->GMAX=source->impulsive->powerlawcutoffelectron->EMAX/electronmassTeV +1.0;
     printf("At the beginning GMIN=%.3e and GMAX=%.3e\n", source->impulsive->powerlawcutoffelectron->GMIN, source->impulsive->powerlawcutoffelectron->GMAX);
    }

    else if(strcmp(source->impulsive->spectratypeelectron,"brokenpowerlawcutoff")==0){
      source->impulsive->brokenpowerlawcutoffelectron=malloc(sizeof(struct brokenpowerlawcutoff)); 
      source->impulsive->brokenpowerlawcutoffelectron=getspectrumbrokenpowerlawcutoff(infile,electronspectrum);
      source->impulsive->brokenpowerlawcutoffelectron->GMIN=source->impulsive->brokenpowerlawcutoffelectron->EMIN/electronmassTeV +1.0;
      source->impulsive->brokenpowerlawcutoffelectron->GMAX=source->impulsive->brokenpowerlawcutoffelectron->EMAX/electronmassTeV +1.0;
    }
  

 /*  proton spectra inut type */
    if (strcmp(source->impulsive->spectratypeproton,"powerlaw")==0){
      source->impulsive->powerlawproton=malloc(sizeof(struct powerlaw));
      source->impulsive->powerlawproton=getspectrumpowerlaw(infile,protonspectrum);
      source->impulsive->powerlawproton->GMIN=source->impulsive->powerlawproton->EMIN/protonmassTeV +1.0;
      source->impulsive->powerlawproton->GMAX=source->impulsive->powerlawproton->EMAX/protonmassTeV +1.0;
    }
    
    else if(strcmp(source->impulsive->spectratypeproton,"brokenpowerlaw")==0){
      source->impulsive->brokenpowerlawproton=malloc(sizeof(struct brokenpowerlaw));
      source->impulsive->brokenpowerlawproton=getspectrumbrokenpowerlaw(infile,protonspectrum);
      source->impulsive->brokenpowerlawproton->GMIN=source->impulsive->brokenpowerlawproton->EMIN/protonmassTeV +1.0;
      source->impulsive->brokenpowerlawproton->GMAX=source->impulsive->brokenpowerlawproton->EMAX/protonmassTeV +1.0;
      
    }
    
    else if(strcmp(source->impulsive->spectratypeproton,"powerlawcutoff")==0){
     source->impulsive->powerlawcutoffproton=malloc(sizeof(struct powerlawcutoff)); 
     source->impulsive->powerlawcutoffproton=getspectrumpowerlawcutoff(infile,protonspectrum);
     source->impulsive->powerlawcutoffproton->GMIN=source->impulsive->powerlawcutoffproton->EMIN/protonmassTeV +1.0;
     source->impulsive->powerlawcutoffproton->GMAX=source->impulsive->powerlawcutoffproton->EMAX/protonmassTeV +1.0;


}

    else if(strcmp(source->impulsive->spectratypeproton,"brokenpowerlawcutoff")==0){
      source->impulsive->brokenpowerlawcutoffproton=malloc(sizeof(struct brokenpowerlawcutoff)); 
      source->impulsive->brokenpowerlawcutoffproton=getspectrumbrokenpowerlawcutoff(infile,protonspectrum);
      source->impulsive->brokenpowerlawcutoffproton->GMIN=source->impulsive->brokenpowerlawcutoffproton->EMIN/protonmassTeV +1;
      source->impulsive->brokenpowerlawcutoffproton->GMAX=source->impulsive->brokenpowerlawcutoffproton->EMAX/protonmassTeV +1;
    }
  




}

  else  if (strcmp(source->sourcetype,"continuous")==0){
    source->continuous=malloc(sizeof(struct continuous));
    sprintf(command,"grep 'Energyrate' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'evolution' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'escape' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'Efractionproton' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'distance' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'size' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'density' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'Bfield' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'time' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'TempIR' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'EnIR' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 

    sprintf(command,"grep 'Tempopt' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'Enopt' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'spectratypeelectron' %s | awk -F '=' '{print $2}' >>tmp.txt",electronspectrum);
    system(command); 
   
   
   
    sprintf(command,"grep 'spectratypeproton' %s | awk -F '=' '{print $2}' >>tmp.txt",protonspectrum);
     system(command);
   
 
    sprintf(command,"grep 'pulsar' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);

    fscanf(infile,"%lf",&(source->continuous->L0));
    fscanf(infile,"%s",source->continuous->evolution);
    fscanf(infile,"%s",source->continuous->escape);
    fscanf(infile,"%lf",&(source->continuous->Efractionproton));
    fscanf(infile,"%lf",&(source->continuous->distance));
    fscanf(infile,"%lf",&(source->continuous->size));
    fscanf(infile,"%lf",&(source->continuous->density));
    fscanf(infile,"%lf",&(source->continuous->Bfield));
    fscanf(infile,"%lf",&(source->continuous->time));
    source->continuous->time *=365.24*3600*24; //convert year to seconds
    fscanf(infile,"%lf",&(source->continuous->TempIR));
    fscanf(infile,"%lf",&(source->continuous->EnergyIR));
    fscanf(infile,"%lf",&(source->continuous->Tempopt));
    fscanf(infile,"%lf",&(source->continuous->Energyopt));
    fscanf(infile,"%s",source->continuous->spectratypeelectron);
    fscanf(infile,"%s",source->continuous->spectratypeproton);
    fscanf(infile,"%s",source->continuous->pulsar);
    printf("%s %s\n",source->continuous->spectratypeelectron,source->continuous->spectratypeproton);
    

    /*electron input spectrum*/
     if (strcmp(source->continuous->spectratypeelectron,"powerlaw")==0){
      source->continuous->powerlawelectron=malloc(sizeof(struct powerlaw));
      source->continuous->powerlawelectron=getspectrumpowerlaw(infile,electronspectrum);
      source->continuous->powerlawelectron->GMIN=source->continuous->powerlawelectron->EMIN/electronmassTeV +1;
      source->continuous->powerlawelectron->GMAX=source->continuous->powerlawelectron->EMAX/electronmassTeV +1;
    }
    
    else if(strcmp(source->continuous->spectratypeelectron,"brokenpowerlaw")==0){
      source->continuous->brokenpowerlawelectron=malloc(sizeof(struct brokenpowerlaw));
      source->continuous->brokenpowerlawelectron=getspectrumbrokenpowerlaw(infile,electronspectrum);
      source->continuous->brokenpowerlawelectron->GMIN=source->continuous->brokenpowerlawelectron->EMIN/electronmassTeV +1;
      source->continuous->brokenpowerlawelectron->GMAX=source->continuous->brokenpowerlawelectron->EMAX/electronmassTeV +1;
    }
    
    else if(strcmp(source->continuous->spectratypeelectron,"powerlawcutoff")==0){
     source->continuous->powerlawcutoffelectron=malloc(sizeof(struct powerlawcutoff)); 
     source->continuous->powerlawcutoffelectron=getspectrumpowerlawcutoff(infile,electronspectrum);
     source->continuous->powerlawcutoffelectron->GMIN=source->continuous->powerlawcutoffelectron->EMIN/electronmassTeV +1;
     source->continuous->powerlawcutoffelectron->GMAX=source->continuous->powerlawcutoffelectron->EMAX/electronmassTeV +1;
     printf("Ola que tal , EMAX=%.3e\n",source->continuous->powerlawcutoffelectron->EMAX);
    }

    else if(strcmp(source->continuous->spectratypeelectron,"brokenpowerlawcutoff")==0){
      source->continuous->brokenpowerlawcutoffelectron=malloc(sizeof(struct brokenpowerlawcutoff)); 
      source->continuous->brokenpowerlawcutoffelectron=getspectrumbrokenpowerlawcutoff(infile,electronspectrum);

      source->continuous->brokenpowerlawcutoffelectron->GMIN=source->continuous->brokenpowerlawcutoffelectron->EMIN/electronmassTeV +1;
      source->continuous->brokenpowerlawcutoffelectron->GMAX=source->continuous->brokenpowerlawcutoffelectron->EMAX/electronmassTeV +1;
      printf("Ok at the beg GMIN=%.3e and GMAX=%.3e \n",source->continuous->brokenpowerlawcutoffelectron->GMIN,source->continuous->brokenpowerlawcutoffelectron->GMAX);
    }

     /*electron input spectrum*/
     if (strcmp(source->continuous->spectratypeproton,"powerlaw")==0){
      
      source->continuous->powerlawproton=malloc(sizeof(struct powerlaw));
      source->continuous->powerlawproton=getspectrumpowerlaw(infile,protonspectrum);
      source->continuous->powerlawproton->GMIN=source->continuous->powerlawproton->EMIN/protonmassTeV +1;
      source->continuous->powerlawproton->GMAX=source->continuous->powerlawproton->EMAX/protonmassTeV +1;
    }
    
    else if(strcmp(source->continuous->spectratypeproton,"brokenpowerlaw")==0){
      source->continuous->brokenpowerlawproton=malloc(sizeof(struct brokenpowerlaw));
      source->continuous->brokenpowerlawproton=getspectrumbrokenpowerlaw(infile,protonspectrum);
      source->continuous->brokenpowerlawproton->GMIN=source->continuous->brokenpowerlawproton->EMIN/protonmassTeV +1;
      source->continuous->brokenpowerlawproton->GMAX=source->continuous->brokenpowerlawproton->EMAX/protonmassTeV +1;
    }
    
    else if(strcmp(source->continuous->spectratypeproton,"powerlawcutoff")==0){
     source->continuous->powerlawcutoffproton=malloc(sizeof(struct powerlawcutoff)); 
     source->continuous->powerlawcutoffproton=getspectrumpowerlawcutoff(infile,protonspectrum);
     source->continuous->powerlawcutoffproton->GMIN=source->continuous->powerlawcutoffproton->EMIN/protonmassTeV +1;
     source->continuous->powerlawcutoffproton->GMAX=source->continuous->powerlawcutoffproton->EMAX/protonmassTeV +1;
   
    }

    else if(strcmp(source->continuous->spectratypeproton,"brokenpowerlawcutoff")==0){
      source->continuous->brokenpowerlawcutoffproton=malloc(sizeof(struct brokenpowerlawcutoff)); 
      source->continuous->brokenpowerlawcutoffproton=getspectrumbrokenpowerlawcutoff(infile,protonspectrum);
      source->continuous->brokenpowerlawcutoffproton->GMIN=source->continuous->brokenpowerlawcutoffproton->EMIN/protonmassTeV +1;
      source->continuous->brokenpowerlawcutoffproton->GMAX=source->continuous->brokenpowerlawcutoffproton->EMAX/protonmassTeV +1;

    }

     /* Now let 's input the pulsar parameters if necessary*/
     if (strcmp(source->continuous->pulsar,"yes")==0){
       source->continuous->pulsartype=malloc(sizeof(struct pulsartype));
       source->continuous->pulsartype=getpulsarparameter(infile,pulsarname);
     }
     
     
     

  }

  else if (strcmp(source->sourcetype,"dirac")==0){
     printf("We are in the dirac option!\n");
   
    source->dirac=malloc(sizeof(struct dirac));
   
    sprintf(command,"grep 'Energy' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'EDirac' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'Bfield' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'density' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'distance' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
     sprintf(command,"grep 'size' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'TempIR' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'EnIR' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'Tempopt' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'Enopt' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    fscanf(infile,"%lf",&(source->dirac->E0));
    fscanf(infile,"%lf",&(source->dirac->Edirac));
    fscanf(infile,"%lf",&(source->dirac->Bfield));
    fscanf(infile,"%lf",&(source->dirac->density));
    fscanf(infile,"%lf",&(source->dirac->distance));
    fscanf(infile,"%lf",&(source->dirac->size));
    fscanf(infile,"%lf",&(source->dirac->TempIR));
    fscanf(infile,"%lf",&(source->dirac->EnergyIR));
    fscanf(infile,"%lf",&(source->dirac->Tempopt));
    fscanf(infile,"%lf",&(source->dirac->Energyopt));
    printf("E0=%3e Edirac=%.3e Bfield=%.2f density=%.2f distance=%.1f size=%.2f TempIr=%.2f EnergyIR=%.2f\n",source->dirac->E0,source->dirac->Edirac,source->dirac->Bfield,source->dirac->density, source->dirac->distance,source->dirac->size,source->dirac->TempIR,source->dirac->EnergyIR  );
  }


  fclose(infile);
}

struct powerlaw *getspectrumpowerlaw(FILE *infile, char filename[200]){
  /*fetch the value*/
  struct powerlaw *tmpspectrum;
  char command[200];
  tmpspectrum=malloc(sizeof(struct powerlaw));
    sprintf(command,"grep 'spectrum1' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
   sprintf(command,"grep 'EMIN' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
  sprintf(command,"grep 'EMAX' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
  system(command);
  
   
  fscanf(infile,"%lf",&(tmpspectrum->inp1));
    fscanf(infile,"%lf",&(tmpspectrum->EMIN));
    fscanf(infile,"%lf",&(tmpspectrum->EMAX));
    printf("EMIN=%.3e EMAX=%.3e\n",tmpspectrum->EMIN,tmpspectrum->EMAX);
    return tmpspectrum;

}

struct brokenpowerlaw *getspectrumbrokenpowerlaw(FILE *infile, char filename[200]){
/*fetch the value*/
  struct brokenpowerlaw *tmpspectrum;
  char command[200];
  tmpspectrum=malloc(sizeof(struct brokenpowerlaw));
   sprintf(command,"grep 'spectrum1' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
  system(command);
  sprintf(command,"grep 'spectrum2' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
  system(command);
   sprintf(command,"grep 'EMIN' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
  sprintf(command,"grep 'EMAX' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);  
  sprintf(command,"grep 'Ebreak' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
  system(command);
  fscanf(infile,"%lf",&(tmpspectrum->inp1));
  fscanf(infile,"%lf",&(tmpspectrum->inp2));
  fscanf(infile,"%lf",&(tmpspectrum->EMIN));
  fscanf(infile,"%lf",&(tmpspectrum->EMAX));
  fscanf(infile,"%lf",&(tmpspectrum->Ebreak));
    return tmpspectrum;

}


struct powerlawcutoff *getspectrumpowerlawcutoff(FILE *infile, char filename[200]){
  /*fetch the value*/
  struct powerlawcutoff *tmpspectrum;
  char command[200];
  tmpspectrum=malloc(sizeof(struct powerlawcutoff));
   sprintf(command,"grep 'spectrum1' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'expspectrum' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
   sprintf(command,"grep 'EMIN' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'EMAX' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);  
    sprintf(command,"grep 'Ecutoff' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    
    
    fscanf(infile,"%lf",&(tmpspectrum->inp1));
    fscanf(infile,"%lf",&(tmpspectrum->expinp));
    fscanf(infile,"%lf",&(tmpspectrum->EMIN));
    fscanf(infile,"%lf",&(tmpspectrum->EMAX));
    fscanf(infile,"%lf",&(tmpspectrum->Ecutoff));
    printf("Did you get there with the right value inp1=%.2f EMAX=%.3e  \n", tmpspectrum->inp1, tmpspectrum->EMAX);
    return tmpspectrum;

}


struct brokenpowerlawcutoff *getspectrumbrokenpowerlawcutoff(FILE *infile, char filename[200]){
  /*fetch the value*/
  struct brokenpowerlawcutoff *tmpspectrum;
  char command[200];
  tmpspectrum=malloc(sizeof(struct brokenpowerlawcutoff));
   sprintf(command,"grep 'spectrum1' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'spectrum2' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'expspectrum' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
   sprintf(command,"grep 'EMIN' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command); 
    sprintf(command,"grep 'EMAX' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);  
    sprintf(command,"grep 'Ecutoff' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
     sprintf(command,"grep 'Ebreak' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    
    
    fscanf(infile,"%lf",&(tmpspectrum->inp1));
    fscanf(infile,"%lf",&(tmpspectrum->inp2));
    fscanf(infile,"%lf",&(tmpspectrum->expinp));
    fscanf(infile,"%lf",&(tmpspectrum->EMIN));
    fscanf(infile,"%lf",&(tmpspectrum->EMAX));
    fscanf(infile,"%lf",&(tmpspectrum->Ecutoff));
    fscanf(infile,"%lf",&(tmpspectrum->Ebreak));
    return tmpspectrum;

}

struct pulsartype *getpulsarparameter(FILE *infile,char filename[200]){
  struct pulsartype *tmppulsar;
  char command[200];
  tmppulsar=malloc(sizeof(struct pulsartype));
  sprintf(command,"grep 'brakingindex' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'period' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'dotPeriod' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'magfraction' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'ejectaMass' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    sprintf(command,"grep 'E_SN' %s | awk -F '=' '{print $2}' >>tmp.txt",filename);
    system(command);
    fscanf(infile,"%lf",&(tmppulsar->brakingindex));
    fscanf(infile,"%lf",&(tmppulsar->period));
    tmppulsar->period *=1E-3; //conversion in s
    fscanf(infile,"%lf",&(tmppulsar->perioddot));
    fscanf(infile,"%lf",&(tmppulsar->fractionmag));
    fscanf(infile,"%lf",&(tmppulsar->Mej));
    tmppulsar->Mej *=solarmass; //conversion in g
    fscanf(infile,"%lf",&(tmppulsar->E_SN));//in 10^51ergs
 
    return tmppulsar;
}

