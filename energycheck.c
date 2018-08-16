#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "../newsedprod/src/constant.h"

int main(){

  int i;
  double step = 5.0;

  int particle = 0;
  double E, dE;
  double integral = 0.0;
  double Wtot;
  double A;
  double mass;
  double Ecut = 0.0;
  double W1TeV = 0.0;
  double alpha1 = 0.0;
  double alpha2 = 0.0;
  double Ebreak = 0.0;
  double Emax = 10000;

  printf("What kind of particle? (1=proton, 2=electron)\n");
  scanf("%d", &particle);
  printf("Enter Wtot, Ecut, Ebreak, alpha1, alpha2 ?\n");
  scanf("%le %le %le %le %le", &Wtot, &Ecut, &Ebreak, &alpha1, &alpha2);
  printf("alpha2 = %.3e\n", alpha2);

  if (particle == 1) mass = protonmassTeV;
  else if (particle == 2) mass = electronmassTeV;
  

  for(i=0;i<250;i++){
    E = 1e-7*pow(10,(double)(i*step)/100.0);
    dE = E * (pow(10,(double)(step)/100.0) - 1.0);
    if((E>mass) & (E<Emax)){
      if(E<Ebreak) integral += pow(E,1) * pow(E,-alpha1) * exp(-(E/Ecut)) * dE;
      else if(E>Ebreak) integral += pow(E,1) * pow(E,-alpha2) * exp(-(E/Ecut)) * dE;
    }
  }

  A = Wtot/integral;

  for(i=0;i<250;i++){
    E = 1e-7*pow(10,(double)(i*step)/100.0);
    dE = E * (pow(10,(step)/100.0) - 1.0);
    if((E>1) & (E<Emax)){
      if(E<Ebreak) W1TeV += A * pow(E,1) * pow(E,-alpha1) * exp(-(E/Ecut)) * dE;
      else if(E>Ebreak) W1TeV += A * pow(E,1) * pow(E,-alpha2) * exp(-(E/Ecut)) * dE;
      } 
  }
  printf("W1TeV = %.3e\n",W1TeV);
  

  return(0);
}
