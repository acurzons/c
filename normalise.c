#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "../newsedprod/src/constant.h"

int main(){

  int i;
  double step = 5.0;

  int type = 0;
  int law = 0;
  int particle = 0;
  int Efactor = 0;
  double E, dE;
  double integral = 0.0;
  double W = 0.0;
  double A;
  double mass = 0.0;
  double threshold_pmassTeV = 1.218e-3;
  double gamma = 0.0;
  double Ecut = 0.0;
  double Fgiven = 0.0;
  double Wgiven = 0.0;
  double gamma2 = 0.0;
  double Ebreak = 0.0;
  double Emax = 0.0;
  double Emin = 0.0;
  double density = 0.0;

  printf("Normalise Particle Energy or Gamma Ray Flux?(1=particle, 2=gammaray flux)\n");
  scanf("%d", &type);

  if(type == 1){
  printf("What kind of particle? (1=proton, 2=electron)\n");
  scanf("%d", &particle);
  printf("What kind of law? (1=cutoff, 2=brokencutoff)\n");
  scanf("%d", &law);
  printf("Given Energy = ?\n");
  scanf("%le", &Wgiven);
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
    if((E>1) & (law == 1)) integral += pow(E,1) * pow(E,-gamma) * exp(-(E/Ecut)) * dE;
    else if((E>1) & (law == 2)){
      if(E<Ebreak) integral += pow(E,1) * pow(E,-gamma) *exp(-(E/Ecut)) * dE;
      else if(E>Ebreak) integral += pow(E,1) * pow(E,-gamma2) * exp(-(E/Ecut)) * dE;
   }
  }

  A = Wgiven/integral;

  for(i=0;i<250;i++){
    E = 1e-7*pow(10,(double)(i*step)/100.0);
    dE = E * (pow(10,(step)/100.0) - 1.0);
    if((E>mass) & (law == 1)) W += A*pow(E,1)*pow(E,-gamma)*exp(-(E/Ecut))*dE;
    else if((E>mass) & (law == 2)){
      if(E<Ebreak) W += A * pow(E,1) * pow(E,-gamma) *exp(-(E/Ecut)) * dE;
      else if(E>Ebreak) W += A * pow(E,1) * pow(E,-gamma2) * exp(-(E/Ecut)) * dE;
      } 
  }
    if(particle == 1) W /= density;
  printf("W = %.3e\n",W);
  }

  else if(type == 2){

    printf("Enter given flux\n");
    scanf("%le", &Fgiven);
    printf("Enter lower energy restriction\n");
    scanf("%le", &Emin);
    printf("Enter upper energy restriction\n");
    scanf("%le", &Emax);
    printf("Enter Gamma\n");
    scanf("%le", &gamma);
    printf("Enter cut off Energy\n");
    scanf("%le", &Ecut);
    printf("Does flux contain a factor of energy? (1 = no, 2 = yes)\n");
    scanf("%d", &Efactor);

    for(i=0; i<2500; i++){
      E = Emin * pow(10,(double)(i*step)/1000.0);
      dE = E *(pow(10,(double)(step)/1000.0) - 1.0);
      if(E<Emax){
	if(Efactor == 1) integral +=  pow(E,-gamma) * exp(-(E/Ecut)) * dE;
	else if(Efactor == 2) integral += pow(E,-gamma+1) * exp(-(E/Ecut)) * dE;
	//      	printf("E = %.3e\n",E);
      }
    }
    A = Fgiven / integral;

    printf("A = %.3e\n", A);
  }

  return(0);
}
