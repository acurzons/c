#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "../newsedprod/src/constant.h"


int main(){

  FILE *fp;
  FILE *params;
  FILE *SED;
  FILE *tmp;
  char file[250], SEDdatafile[250], buffer[250], line[100], command[100], paramsline[100];
  char paramsfile[250] = "rxj1713_regional_hadronic_parameters.txt";
  
  int i, j, N, reg[29];
 
  double E[29][100], flux[29][100], fluxerror[29][100], n[29], Wpm[29], alpha[29], Ec[29], fod;

  double Eint;
  double dEint;
  double step = 5.0;
  double A;
  double integral = 0.0;
  double Wlimits[29];
  double distance = 1000*3.1e18; //cm

  double *ESED[29];
  double *fluxSED[29];
  double fluxsumSED[29];
  double LpSED[29];
  double tSED[29];

  double Emin[29], Emax[29];
  double dE;
  double avedE;
  double aveE;
  double aveflux;
  double avefluxsum[29];
  double fluxsum[29];
  double Lp[29];
  double Wpt[29];
  double aveWpt[29];
  double tpp[29];

  double ratio[29];
  double averatio[29];

  double tppyears[29];
  double t[29];
  double tremaining[29];
  double tpercent[29];

  for(j=0;j<29;j++){
    Wlimits[j] = 0.0;
    for(i=0;i<100;i++){
      E[j][i] = 0.0;
      flux[j][i] = 0.0;
    }
  }

  //Read in hadronic SED parameters 
  params = fopen(paramsfile, "r");
  fgets(buffer, 200, params);
  for(j=0;j<29;j++){
    fgets(paramsline, sizeof(paramsline), params);
    sscanf(paramsline, "%d %lf %lf %lf %lf %lf %lf %lf", &reg[j], &fod, &n[j], &Wpm[j], &fod, &alpha[j], &Ec[j], &fod);
    Wpm[j] *= 1e46;


    //Open a temp file to write how many data points in each SED data file
    sprintf(SEDdatafile, "../newsedprod/OUTPUT/rxj1713_had%d_ppgamma.txt", j+1);
    tmp = fopen("tmp.txt", "r");
    sprintf(command, "cat %s | wc -l > tmp.txt", SEDdatafile);
    system(command);
    fscanf(tmp, "%d", &N);
    fclose(tmp);
  
    //Allocate memory
    ESED[j] = malloc(sizeof(double) * N);
    fluxSED[j] = malloc(sizeof(double) * N);

    //Read in SED data
    SED = fopen(SEDdatafile, "r");
    i = 0;
    while(fgets(line, sizeof(line), SED)) {
      sscanf(line, "%lf %lf %lf %lf %lf %lf", &ESED[j][i], &fod, &fod, &fod, &fod, &fluxSED[j][i]);
      i++;
    }
  
    //Compute total flux - area under SED curve (convert everything to erg *1.6)
    for(i=0;i<N;i++){
      dE = (ESED[j][i+1] - ESED[j][i]) * 1.6;
      fluxsumSED[j] += fluxSED[j][i] * dE / (ESED[j][i] * 1.6);
    }
    
    //Compute age of SNR in years
    LpSED[j] = 4 * PI * pow(distance,2) * fluxsumSED[j];
    tSED[j] = Wpm[j] / LpSED[j] / (60*60*24*365);     
    
    sprintf(file, "../newsedprod/OUTPUT/authorsplots/hess_spectra/region_%d.spec", j+1);
    fp = fopen(file, "r");
    fgets(buffer, 100, fp);
    fgets(buffer, 100, fp);
    fgets(buffer, 100, fp);

    i = 0;
    while(fgets(line, sizeof(line), fp)) {
      sscanf(line, "%lf %lf %lf", &E[j][i], &flux[j][i], &fluxerror[j][i]);
      //      if (j == 0) printf("E = %.3lf, i = %d\n",E[j][i], i);
      i++;
    } 

    Emin[j] = E[j][0];
    Emax[j] = E[j][0];
    i = 0;
    Lp[j] = 0.0;
    while(E[j][i+1]!=0){       
      if(E[j][i+1]<Emin[j]) Emin[j] = E[j][i+1];
      if(E[j][i+1]>Emax[j]) Emax[j] = E[j][i+1];
      dE = E[j][i+1] - E[j][i];
      dE *= 1.6; //TeV to erg
      fluxsum[j] += flux[j][i] * dE / (E[j][i]*1.6); 
      i++;
    }

    i = 0;
    while(E[j][i+1]!=0){
      aveE = (E[j][i] + E[j][i+1]) * 1.6 / 2.0;
      dE = aveE - E[j][i] * 1.6;
      avedE = E[j][i+1]*1.6 - aveE; 
      aveflux = (flux[j][i] + flux[j][i+1]) / 2.0;
      avefluxsum[j] += flux[j][i] * dE / (E[j][i]*1.6) +  aveflux * avedE / aveE;
      i++;
    }


    Lp[j] = 4 * PI * pow(distance,2) * fluxsum[j];
    tpp[j] = 5.3e7 / n[j] * 365 * 24 * 60 * 60;
    Wpt[j] = Lp[j] * tpp[j];

    aveWpt[j] = 4 * PI * pow(distance,2) * avefluxsum[j] * tpp[j];

    for(i=0;i<250;i++){
      Eint = 1e-7*pow(10,(double)(i*step)/100.0);
      dEint = Eint * (pow(10,(double)(step)/100.0) - 1.0);
      if(Eint>1.2e-3) integral += pow(Eint,1) * pow(Eint,-alpha[j]) * exp(-(Eint/Ec[j])) * dEint;
      }
    
    A = Wpm[j]/integral;
    
    for(i=0;i<250;i++){
      Eint = 1e-7*pow(10,(double)(i*step)/100.0);
      dEint = Eint * (pow(10,(step)/100.0) - 1.0);
      if((Eint>Emin[j]) & (Eint<Emax[j])) Wlimits[j] += A*pow(Eint,1)*pow(Eint,-alpha[j])*exp(-(Eint/Ec[j]))*dEint;
      }

    ratio[j] = Wpt[j] / Wlimits[j];
    averatio[j] = aveWpt[j] / Wlimits[j];

    tppyears[j] = tpp[j]/(60*60*24*365);
    t[j] = Wlimits[j] / Lp[j] / (60*60*24*365);
    tremaining[j] = tppyears[j] - t[j];
    tpercent[j] = tppyears[j] / tremaining[j];

    printf("R%d, Wpt = %.2e, Wpm = %.2e, ratio = %.1lf, tSED = %.2e, tF = %.2e, tpp = %.2e\n",reg[j], Wpt[j], Wlimits[j], ratio[j], tSED[j], t[j], tppyears[j]);



  }
  
  fclose(fp);
  return(0);

}
