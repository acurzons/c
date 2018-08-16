#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "../newsedprod/src/constant.h"

int main(){

  int i;
  double step = 5.0;

  int law = 0;
  int particle = 0;
  double E, dE;
  double integral = 0.0;
  double W = 0.0;
  double A;
  double mass = 0.0;
  double threshold_pmassTeV = 1.218e-3;
  double gamma = 0.0;
  double Ecut = 0.0;
  double gamma2 = 0.0;
  double Ebreak = 0.0;
  double Emax = 0.0;
  double Emin = 0.0;
  double density = 0.0;
  double Wlimits = 0.0;


  printf("What kind of particle? (1=proton, 2=electron)\n");
  scanf("%d", &particle);
  printf("What kind of law? (1=cutoff, 2=brokencutoff)\n");
  scanf("%d", &law);
  printf("Total Injection Energy = ?\n");
  scanf("%le", &W);
  printf("Emin = ?\n");
  scanf("%le", &Emin);
  printf("Emax = ?\n");
  scanf("%le", &Emax);

  //  Wgiven *= ergtoTeV;
  printf("Ecut = ?\n");
  scanf("%le", &Ecut);
  printf("gamma = ?\n");
  scanf("%le", &gamma);
  if(law == 2){
    printf("gamma2 = ?\n");
    scanf("%le", &gamma2);
    printf("Ebreak = ?\n");
    scanf("%le", &Ebreak);
  }

  if (particle == 1){
    mass = threshold_pmassTeV;
    //    mass = protonmassTeV;
    printf("What is the denisty?\n");
    scanf("%le", &density);
  }
  else if (particle == 2) mass = electronmassTeV;
  
    printf("mass = %.3e\n", mass);

  for(i=0;i<250;i++){
    E = 1e-7*pow(10,(double)(i*step)/100.0);
    dE = E * (pow(10,(double)(step)/100.0) - 1.0);
    //    printf("E = %.3e, dE = %.3e\n", E, dE);
    if((E>mass) & (law == 1)) integral += pow(E,1) * pow(E,-gamma) * exp(-(E/Ecut)) * dE;
    else if((E>mass) & (law == 2)){
      if(E<Ebreak) integral += pow(E,1) * pow(E,-gamma) *exp(-(E/Ecut)) * dE;
      else if(E>Ebreak) integral += pow(E,1) * pow(E,-gamma2) * exp(-(E/Ecut)) * dE;
    }
  }

  A = W/integral;
  printf("A = %.3e\n",A);
  for(i=0;i<250;i++){
    E = 1e-7*pow(10,(double)(i*step)/100.0);
    dE = E * (pow(10,(step)/100.0) - 1.0);
    if((E>Emin) & (E<Emax) & (law == 1)) Wlimits += A*pow(E,1)*pow(E,-gamma)*exp(-(E/Ecut))*dE;
    else if((E>Emin) & (E<Emax) & (law == 2)){
      if(E<Ebreak) Wlimits += A * pow(E,1) * pow(E,-gamma) *exp(-(E/Ecut)) * dE;
      else if(E>Ebreak) Wlimits += A * pow(E,1) * pow(E,-gamma2) * exp(-(E/Ecut)) * dE;
    } 
  }
  printf("W in your energy range = %.2e\n", Wlimits);


  return(0);
}
