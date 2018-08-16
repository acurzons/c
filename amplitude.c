#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "../newsedprod/src/constant.h"


int main(){

  FILE *fp;
  FILE *params;
  char file[250];
  char paramsfile[250] = "rxj1713_regional_hadronic_parameters.txt";
  char buffer[250];
  char line[100];
  char paramsline[100];

  int i;
  int j;
 
  double E[29][100];
  double flux[29][100];
  double fluxerror[29][100];
  int reg[29];
  double n[29];
  double Wpm[29];
  double alpha[29];
  double Ec[29];
  double fod;

  double Eint;
  double dEint;
  double step = 5.0;
  double A[29];
  double integral = 0.0;
  double distance = 1000*3.1e18; //cm
 
  params = fopen(paramsfile, "r");
  fgets(buffer, 200, params);

  for(j=0;j<29;j++){

    fgets(paramsline, sizeof(paramsline), params);
    sscanf(paramsline, "%d %lf %lf %lf %lf %lf %lf %lf", &reg[j], &fod, &n[j], &Wpm[j], &fod, &alpha[j], &Ec[j], &fod);
    Wpm[j] *= 1e46;
    
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

    i = 0;

    for(i=0;i<250;i++){
      Eint = 1e-7*pow(10,(double)(i*step)/100.0);
      dEint = Eint * (pow(10,(double)(step)/100.0) - 1.0);
      if(Eint>1.2e-3) integral += pow(Eint,1) * pow(Eint,-alpha[j]) * exp(-(Eint/Ec[j])) * dEint;
      }
    
    A[j] = (Wpm[j]/1.6)/integral;
    A[j] /= 1.0e12; //convert to eV

    printf("Reg.%d, A = %.3e\n", reg[j], A[j]);
  }
  fclose(fp);
  return(0);

}