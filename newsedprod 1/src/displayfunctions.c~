#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structure.h"
#include "constant.h"
#include "globalparameters.h"
#include "displayfunctionslib.h"
#include "pulsarmisclib.h"


void makeoutput(char prefix[50]){

  int i;
  double size,distance;
  if(source->continuous!=NULL){
    if (source->continuous->pulsartype==NULL) size=source->continuous->size*parsectocm;
    else {
      printf("source->time =%.3e\n",source->continuous->time);
      size=getsize(source->continuous->time);
      printf("fixglitch\n");
      printf("size3=%.3e\n",size/parsectocm);
    }
    distance=source->continuous->distance*parsectocm;
  }
  else if (source->impulsive!=NULL){
    size=source->impulsive->size*parsectocm;
    distance=source->impulsive->distance*parsectocm;
  }
  else if (source->dirac!=NULL) {
    size=source->dirac->size*parsectocm;
    distance=source->dirac->distance*parsectocm;
  }

  printf("size=%.3e\n",size);

  FILE *syncprimfile;
  FILE *ICprimfile;
  FILE *ICRprimfile;
  FILE *ICOprimfile;
  FILE *bremssprimfile;
  FILE *syncsecondaryfile;
  FILE *secondaryelectronfile;
  FILE *bremsssecondaryfile;
  FILE *SSCfile;
  FILE *ppgammafile;
  FILE *total;
  FILE *electron_dist;
  FILE *proton_dist;

  char syncprimname[100];
  char ICprimname[100];
  char ICRprimname[100];
  char ICOprimname[100];
  char bremssprimname[100];
  char syncsecondaryname[100];
  char bremsssecondaryname[100];
  char secondaryelectronname[100];
  char SSCname[100];
  char ppgammaname[100];
  char totalname[100];
  

  char electronname[100];
  char protonname[100];
  sprintf(syncprimname,"%sOUTPUT/%s_syncprimary.txt",getenv("SEDPRODFOLDER"),prefix);
  sprintf(ICprimname,"%sOUTPUT/%s_ICprimary.txt",getenv("SEDPRODFOLDER"),prefix);
  sprintf(ICRprimname,"%sOUTPUT/%s_ICRprimary.txt",getenv("SEDPRODFOLDER"),prefix);
   sprintf(ICOprimname,"%sOUTPUT/%s_ICOprimary.txt",getenv("SEDPRODFOLDER"),prefix);
  sprintf(bremssprimname,"%sOUTPUT/%s_bremssprimary.txt",getenv("SEDPRODFOLDER"),prefix);
  sprintf(syncsecondaryname,"%sOUTPUT/%s_syncsecondary.txt",getenv("SEDPRODFOLDER"),prefix);
  sprintf(bremsssecondaryname,"%sOUTPUT/%s_bremsssecondary.txt",getenv("SEDPRODFOLDER"),prefix);
   sprintf(SSCname,"%sOUTPUT/%s_SSC.txt",getenv("SEDPRODFOLDER"),prefix);
   sprintf(ppgammaname,"%sOUTPUT/%s_ppgamma.txt",getenv("SEDPRODFOLDER"),prefix);
    sprintf(secondaryelectronname,"%sOUTPUT/%s_secondaryelec.txt",getenv("SEDPRODFOLDER"),prefix);
   sprintf(totalname,"%sOUTPUT/%s_TOTAL.txt",getenv("SEDPRODFOLDER"),prefix);
   sprintf(protonname,"%sOUTPUT/%s_proton.txt",getenv("SEDPRODFOLDER"),prefix);
   sprintf(electronname,"%sOUTPUT/%s_electron.txt",getenv("SEDPRODFOLDER"),prefix);

   if(source->dirac==NULL){
     /*creation of the different file*/
    
     if((syncprimfile=fopen(syncprimname,"w+"))==NULL){
       printf("unable to open %s",syncprimname);
       exit(1);
     }

     if((ICprimfile=fopen(ICprimname,"w+"))==NULL){
       printf("unable to open %s",ICprimname);
       exit(1);
     }
     
     if((ICRprimfile=fopen(ICRprimname,"w+"))==NULL){
       printf("unable to open %s",ICRprimname);
       exit(1);
     }

     if((ICOprimfile=fopen(ICOprimname,"w+"))==NULL){
       printf("unable to open %s",ICRprimname);
       exit(1);
     }

     if((bremssprimfile=fopen(bremssprimname,"w+"))==NULL){
       printf("unable to open %s",bremssprimname);
       exit(1);
     }

     if((syncsecondaryfile=fopen(syncsecondaryname,"w+"))==NULL){
       printf("unable to open %s",syncsecondaryname);
       exit(1);
     }

      if((bremsssecondaryfile=fopen(bremsssecondaryname,"w+"))==NULL){
       printf("unable to open %s",bremsssecondaryname);
       exit(1);
     }


       if((SSCfile=fopen(SSCname,"w+"))==NULL){
       printf("unable to open %s",SSCname);
       exit(1);
     }

       if((ppgammafile=fopen(ppgammaname,"w+"))==NULL){
       printf("unable to open %s",ppgammaname);
       exit(1);
     }

        if((ppgammafile=fopen(ppgammaname,"w+"))==NULL){
       printf("unable to open %s",ppgammaname);
       exit(1);
     }
     
	if((secondaryelectronfile=fopen(secondaryelectronname,"w+"))==NULL){
       printf("unable to open %s",ppgammaname);
       exit(1);
     }
       if((total=fopen(totalname,"w+"))==NULL){
       printf("unable to open %s",totalname);
       exit(1);
     }
        if((electron_dist=fopen(electronname,"w+"))==NULL){
       printf("unable to open %s",electronname);
       exit(1);
     }
	
	if((proton_dist=fopen(protonname,"w+"))==NULL){
       printf("unable to open %s",protonname);
       exit(1);
	}
	totalphoton();
       /*ALL THE PHOTON FILE ARE SUPPOSED TO HAVE THE SAME SIZE*/
	/*the format will be as follow :
	  -Energy in TeV
	  -photon flux in ph s-1 TeV-1
	  -Energy flux distribution in erg s-1
	  -photon energy density in erg cm-3 s-1
	  -Energy flux distribution received in erg cm-2 s-1 str-1
	  
	*/
       for(i=0;i<syncgamma->size;i++){
	 fprintf(syncprimfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",syncgamma->array[0][i],syncgamma->array[1][i],1.0/ergtoTeV*pow(syncgamma->array[0][i],2)*syncgamma->array[1][i],1.0/ergtoTeV*pow(syncgamma->array[0][i],2)*syncgamma->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(4.0*PI)/ergtoTeV*pow(syncgamma->array[0][i],2)*syncgamma->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(syncgamma->array[0][i],2)*syncgamma->array[1][i]/(4.0*PI*pow(distance,2.0)));

	 fprintf(ICprimfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",ICgamma->array[0][i],ICgamma->array[1][i],1.0/ergtoTeV*pow(ICgamma->array[0][i],2)*ICgamma->array[1][i],1.0/ergtoTeV*pow(ICgamma->array[0][i],2)*ICgamma->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(2.0*PI)/ergtoTeV*pow(ICgamma->array[0][i],2)*ICgamma->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(ICgamma->array[0][i],2)*ICgamma->array[1][i]/(4.0*PI*pow(distance,2.0)));
	 
	 fprintf(ICRprimfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",ICRgamma->array[0][i],ICRgamma->array[1][i],1.0/ergtoTeV*pow(ICRgamma->array[0][i],2)*ICRgamma->array[1][i],1.0/ergtoTeV*pow(ICRgamma->array[0][i],2)*ICRgamma->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(2.0*PI)/ergtoTeV*pow(ICRgamma->array[0][i],2)*ICRgamma->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(ICRgamma->array[0][i],2)*ICRgamma->array[1][i]/(4.0*PI*pow(distance,2.0)));

 	fprintf(ICOprimfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",OPTgamma->array[0][i],OPTgamma->array[1][i],1.0/ergtoTeV*pow(OPTgamma->array[0][i],2)*OPTgamma->array[1][i],1.0/ergtoTeV*pow(OPTgamma->array[0][i],2)*OPTgamma->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(2.0*PI)/ergtoTeV*pow(OPTgamma->array[0][i],2)*OPTgamma->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(OPTgamma->array[0][i],2)*OPTgamma->array[1][i]/(4.0*PI*pow(distance,2.0)));
	 
	 fprintf(bremssprimfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",bremssgamma->array[0][i],bremssgamma->array[1][i],1.0/ergtoTeV*pow(bremssgamma->array[0][i],2)*bremssgamma->array[1][i],1.0/ergtoTeV*pow(bremssgamma->array[0][i],2)*bremssgamma->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(4.0*PI)/ergtoTeV*pow(bremssgamma->array[0][i],2)*bremssgamma->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(bremssgamma->array[0][i],2)*bremssgamma->array[1][i]/(4.0*PI*pow(distance,2.0)));


	  fprintf(syncsecondaryfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",secondaryelecgammasync->array[0][i],secondaryelecgammasync->array[1][i],1.0/ergtoTeV*pow(secondaryelecgammasync->array[0][i],2)*secondaryelecgammasync->array[1][i],1.0/ergtoTeV*pow(secondaryelecgammasync->array[0][i],2)*secondaryelecgammasync->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(4.0*PI)/ergtoTeV*pow(secondaryelecgammasync->array[0][i],2)*secondaryelecgammasync->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(secondaryelecgammasync->array[0][i],2)*secondaryelecgammasync->array[1][i]/(4.0*PI*pow(distance,2.0)));
	 
	  fprintf(bremsssecondaryfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",secondaryelecgammabremss->array[0][i],secondaryelecgammabremss->array[1][i],1.0/ergtoTeV*pow(secondaryelecgammabremss->array[0][i],2)*secondaryelecgammabremss->array[1][i],1.0/ergtoTeV*pow(secondaryelecgammabremss->array[0][i],2)*secondaryelecgammabremss->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(4.0*PI)/ergtoTeV*pow(secondaryelecgammabremss->array[0][i],2)*secondaryelecgammabremss->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(secondaryelecgammabremss->array[0][i],2)*secondaryelecgammabremss->array[1][i]/(4.0*PI*pow(distance,2.0)));

	  if (SSCoption!=0){
	    fprintf(SSCfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",SSCgamma->array[0][i],SSCgamma->array[1][i],1.0/ergtoTeV*pow(SSCgamma->array[0][i],2)*SSCgamma->array[1][i],1.0/ergtoTeV*pow(SSCgamma->array[0][i],2)*SSCgamma->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(4.0*PI)/ergtoTeV*pow(SSCgamma->array[0][i],2)*SSCgamma->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(SSCgamma->array[0][i],2)*SSCgamma->array[1][i]/(4.0*PI*pow(distance,2.0)));
	  }
	  fprintf(ppgammafile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",gammaraypp->array[0][i],gammaraypp->array[1][i],1.0/ergtoTeV*pow(gammaraypp->array[0][i],2)*gammaraypp->array[1][i],1.0/ergtoTeV*pow(gammaraypp->array[0][i],2)*gammaraypp->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(4.0*PI)/ergtoTeV*pow(gammaraypp->array[0][i],2)*gammaraypp->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(gammaraypp->array[0][i],2)*gammaraypp->array[1][i]/(4.0*PI*pow(distance,2.0)));


	  fprintf(secondaryelectronfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",secondaryelecpp->array[0][i],secondaryelecpp->array[1][i],1.0/ergtoTeV*pow(secondaryelecpp->array[0][i],2)*secondaryelecpp->array[1][i],1.0/ergtoTeV*pow(secondaryelecpp->array[0][i],2)*secondaryelecpp->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(4.0*PI)/ergtoTeV*pow(secondaryelecpp->array[0][i],2)*secondaryelecpp->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(secondaryelecpp->array[0][i],2)*secondaryelecpp->array[1][i]/(4.0*PI*pow(distance,2.0)));

	  fprintf(total,"%.3e %.3e %.3e %.3e %.3e %.3e\n",totalgamma->array[0][i],totalgamma->array[1][i],1.0/ergtoTeV*pow(totalgamma->array[0][i],2)*totalgamma->array[1][i],1.0/ergtoTeV*pow(totalgamma->array[0][i],2)*totalgamma->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(4.0*PI)/ergtoTeV*pow(totalgamma->array[0][i],2)*totalgamma->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(totalgamma->array[0][i],2)*totalgamma->array[1][i]/(4.0*PI*pow(distance,2.0)));

	   



	   
       }
      
       /*electron distribution*/
      

       for(i=0;i<electronarraysize;i++){
	 fprintf(electron_dist,"%.3e %.3e %.3e\n",electrondistribution[0][i]*electronmassTeV,electrondistribution[1][i],pow(electrondistribution[0][i]*electronmassTeV,2)*electrondistribution[1][i]);
       }

      
       for(i=0;i<protonarraysize;i++){
	 fprintf(proton_dist,"%.3e %.3e %.3e\n",protondistribution[0][i]*protonmassTeV,protondistribution[1][i],pow(protondistribution[0][i]*protonmassTeV,2)*protondistribution[1][i]);

       }
       /*closing all the file*/

       fclose(syncprimfile);
       fclose(ICprimfile);
       fclose(ICRprimfile);
       fclose(ICOprimfile);
       fclose(bremssprimfile);
       fclose(syncsecondaryfile);
       fclose(bremsssecondaryfile);
       fclose(SSCfile);
       fclose(ppgammafile);
       fclose(total);
       fclose(electron_dist);
       fclose(proton_dist);



   }


   else if(source->dirac!=NULL){
    
   if((syncprimfile=fopen(syncprimname,"w+"))==NULL){
       printf("unable to open %s",syncprimname);
       exit(1);
     }

     if((ICprimfile=fopen(ICprimname,"w+"))==NULL){
       printf("unable to open %s",ICprimname);
       exit(1);
     }
     
     if((ICRprimfile=fopen(ICRprimname,"w+"))==NULL){
       printf("unable to open %s",ICRprimname);
       exit(1);
     }
     if((ICOprimfile=fopen(ICOprimname,"w+"))==NULL){
       printf("unable to open %s",ICRprimname);
       exit(1);
     }

     if((bremssprimfile=fopen(bremssprimname,"w+"))==NULL){
       printf("unable to open %s",bremssprimname);
       exit(1);
     }
  
      if((ppgammafile=fopen(ppgammaname,"w+"))==NULL){
       printf("unable to open %s",ppgammaname);
       exit(1);
     }
       if((total=fopen(totalname,"w+"))==NULL){
       printf("unable to open %s",totalname);
       exit(1);
     }
 
      totalphotondirac();
      printf("DIRAC OUTPUT!!!\n");

 for(i=0;i<syncgamma->size;i++){
   fprintf(syncprimfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",syncgamma->array[0][i],syncgamma->array[1][i],ergtoTeV*pow(syncgamma->array[0][i],2)*syncgamma->array[1][i],ergtoTeV*pow(syncgamma->array[0][i],2)*syncgamma->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(4.0*PI)/ergtoTeV*pow(syncgamma->array[0][i],2)*syncgamma->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(syncgamma->array[0][i],2)*syncgamma->array[1][i]/(4.0*PI*pow(distance,2.0)));

   fprintf(ICprimfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",ICgamma->array[0][i],ICgamma->array[1][i],ergtoTeV*pow(ICgamma->array[0][i],2)*ICgamma->array[1][i],ergtoTeV*pow(ICgamma->array[0][i],2)*ICgamma->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(2.0*PI)/ergtoTeV*pow(ICgamma->array[0][i],2)*ICgamma->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(ICgamma->array[0][i],2)*ICgamma->array[1][i]/(4.0*PI*pow(distance,2.0)));
	 
   fprintf(ICRprimfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",ICRgamma->array[0][i],ICRgamma->array[1][i],ergtoTeV*pow(ICRgamma->array[0][i],2)*ICRgamma->array[1][i],ergtoTeV*pow(ICRgamma->array[0][i],2)*ICRgamma->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(2.0*PI)/ergtoTeV*pow(ICRgamma->array[0][i],2)*ICRgamma->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(ICRgamma->array[0][i],2)*ICRgamma->array[1][i]/(4.0*PI*pow(distance,2.0)));

   fprintf(ICOprimfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",OPTgamma->array[0][i],OPTgamma->array[1][i],ergtoTeV*pow(OPTgamma->array[0][i],2)*OPTgamma->array[1][i],ergtoTeV*pow(OPTgamma->array[0][i],2)*OPTgamma->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(2.0*PI)/ergtoTeV*pow(OPTgamma->array[0][i],2)*OPTgamma->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(OPTgamma->array[0][i],2)*OPTgamma->array[1][i]/(4.0*PI*pow(distance,2.0)));
	 
   fprintf(bremssprimfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",bremssgamma->array[0][i],bremssgamma->array[1][i],ergtoTeV*pow(bremssgamma->array[0][i],2)*bremssgamma->array[1][i],ergtoTeV*pow(bremssgamma->array[0][i],2)*bremssgamma->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(4.0*PI)/ergtoTeV*pow(bremssgamma->array[0][i],2)*bremssgamma->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(bremssgamma->array[0][i],2)*bremssgamma->array[1][i]/(4.0*PI*pow(distance,2.0)));

   fprintf(ppgammafile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",gammaraypp->array[0][i],gammaraypp->array[1][i],ergtoTeV*pow(gammaraypp->array[0][i],2)*gammaraypp->array[1][i],ergtoTeV*pow(gammaraypp->array[0][i],2)*gammaraypp->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(4.0*PI)/ergtoTeV*pow(gammaraypp->array[0][i],2)*gammaraypp->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(gammaraypp->array[0][i],2)*gammaraypp->array[1][i]/(4.0*PI*pow(distance,2.0)));

   fprintf(total,"%.3e %.3e %.3e %.3e %.3e %.3e\n",totalgamma->array[0][i],totalgamma->array[1][i],ergtoTeV*pow(totalgamma->array[0][i],2)*totalgamma->array[1][i],ergtoTeV*pow(totalgamma->array[0][i],2)*totalgamma->array[1][i]/(4.0/3.0*PI*pow(size,3.0)),1.0/(4.0*PI)/ergtoTeV*pow(totalgamma->array[0][i],2)*totalgamma->array[1][i]/(4.0*PI*pow(distance,2.0)),1.0/ergtoTeV*pow(totalgamma->array[0][i],2)*totalgamma->array[1][i]/(4.0*PI*pow(distance,2.0)));




 }

 fclose(syncprimfile);
 fclose(ICprimfile);
 fclose(ICRprimfile);
 fclose(ICOprimfile);
 fclose(bremssprimfile);
 fclose(ppgammafile);
 fclose(total);


   }

}



void totalphoton(){
  /*return an array with the total number of photon*/
  int i,j;
  
  totalgamma=malloc(sizeof(struct ARRAYSIZE));
  
  totalgamma->size=syncgamma->size;
  for(i=0;i<2;i++){
    totalgamma->array[i]=malloc(totalgamma->size*sizeof(double));
    for(j=0;j<totalgamma->size;j++){
     totalgamma->array[i][j]=0.0;
    }
  }
  printf("test\n");
  for(i=0;i<totalgamma->size;i++){
    totalgamma->array[0][i]=syncgamma->array[0][i];
    if(SSCoption==1){
    totalgamma->array[1][i]=gammaraypp->array[1][i]+secondaryelecgammasync->array[1][i]+secondaryelecgammabremss->array[1][i]+ICgamma->array[1][i]+OPTgamma->array[1][i]+syncgamma->array[1][i]+bremssgamma->array[1][i]+ICRgamma->array[1][i]+SSCgamma->array[1][i];
    }
    else{
       totalgamma->array[1][i]=gammaraypp->array[1][i]+secondaryelecgammasync->array[1][i]+secondaryelecgammabremss->array[1][i]+ICgamma->array[1][i]+syncgamma->array[1][i]+bremssgamma->array[1][i]+ICRgamma->array[1][i]+OPTgamma->array[1][i];
    }
    /*sum of all gamma SED*/


  }
    
}





void totalphotondirac(){
  /*return an array with the total number of photon*/
  int i,j;
  
  totalgamma=malloc(sizeof(struct ARRAYSIZE));
  
  totalgamma->size=syncgamma->size;
  for(i=0;i<2;i++){
    totalgamma->array[i]=malloc(totalgamma->size*sizeof(double));
    for(j=0;j<totalgamma->size;j++){
     totalgamma->array[i][j]=0.0;
    }
  }
  printf("test\n");
  for(i=0;i<totalgamma->size;i++){
    totalgamma->array[0][i]=syncgamma->array[0][i];
    if(SSCoption==1){
    totalgamma->array[1][i]=gammaraypp->array[1][i]+ICgamma->array[1][i]+syncgamma->array[1][i]+bremssgamma->array[1][i]+ICRgamma->array[1][i]+SSCgamma->array[1][i]+OPTgamma->array[1][i];
    }
    else{
       totalgamma->array[1][i]=gammaraypp->array[1][i]+ICgamma->array[1][i]+syncgamma->array[1][i]+bremssgamma->array[1][i]+ICRgamma->array[1][i]+OPTgamma->array[1][i];
    }
    /*sum of all gamma SED*/


  }
    
}
