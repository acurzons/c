#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constant.h"
#include "structure.h"
#include "globalparameters.h"
#include "pulsarmisclib.h"

double getmagneticfield(double time){
  /*this function will read the pulsarBfield array and indicate its Magnetic field at this epoch*/
  int i=0;
  double B;
  if (magevolution==0){
    B=source->continuous->Bfield; /*already in microgauss*/
    //printf("preBfield\n");
    return B;
  }
  else{
  while(time>pulsarBfield[0][i]){
    //printf("pulsarBfield1=%.3e\n",pulsarBfield[0][i]);
    i++;
  }
  B=pulsarBfield[1][i];
  //printf("Bfieldpick=%.3e\n",B*1E6);
 
  return B*1E6; //MUST RETURN THE MAGNETIC FIELD IN MICROGAUSS
  }



}






double rk4(double(*f)(double, double), double dx, double x, double y)
{
	double	k1 = dx * f(x, y),
		k2 = dx * f(x + dx / 2, y + k1 / 2),
		k3 = dx * f(x + dx / 2, y + k2 / 2),
		k4 = dx * f(x + dx, y + k3);
	return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}
 
double Bfield(double Wb, double time)
{
  
  
  double Bfield=sqrt(6*Wb/(pow(C1*Vej,3.0)*pow(time,18.0/5.0)));
  
  return Bfield;
  
	
}
 

double rate(double t, double y)
{
  double P0=source->continuous->pulsartype->period*pow(1-source->continuous->pulsartype->perioddot*(source->continuous->pulsartype->brakingindex-1.0)*source->continuous->time/(source->continuous->pulsartype->period),1.0/(source->continuous->pulsartype->brakingindex-1));
  double Pdot0=pow(source->continuous->pulsartype->period/P0,source->continuous->pulsartype->brakingindex-2)*source->continuous->pulsartype->perioddot;
  double t0=P0/((source->continuous->pulsartype->brakingindex-1)*Pdot0);
  printf("t0=%.3e, P0=%.3e Pdot=%.3e\n",t0,P0,Pdot0);
  return source->continuous->pulsartype->fractionmag*source->continuous->L0*pow(1.0+t/t0,-(source->continuous->pulsartype->brakingindex+1.0)/(source->continuous->pulsartype->brakingindex-1.0))-6.0/5.0*y/t;
  //return eta*L0-6.0/5.0*y/t;

}

double initialcondition(double time){
 
  return 5.0/11.0*source->continuous->pulsartype->fractionmag*source->continuous->L0*time; 



}
void getBfieldlookuptable()
{
	double *y, x, y2;
       
	double x0 = 5*yeartosec, x1 = source->continuous->time+100*yeartosec, dx =10*yeartosec;
	int i, n = 1 + (x1 - x0)/dx;
	
	C1=pow(6.0/5.0+289/240,-0.2)*pow(source->continuous->L0/source->continuous->pulsartype->E_SN,0.2);
	Vej=sqrt(10.0*source->continuous->pulsartype->E_SN/(3.0*source->continuous->pulsartype->Mej));
	
	printf("C1=%.3e  Vej=%.3e E_SN=%.3e\n",C1,Vej,source->continuous->pulsartype->E_SN);
	

	y = malloc(sizeof(double) * n);
	/*initialization of the look up array*/
	for(i=0;i<2;i++){
	  pulsarBfield[i]=malloc(n*sizeof(double));
	}
	for (y[0] = initialcondition(x0), i = 1; i < n; i++){
	  y[i] = rk4(rate, dx, x0 + dx * (i - 1), y[i-1]);
	  pulsarBfield[0][i]=x0+i*dx;
	  pulsarBfield[1][i]=Bfield(y[i],x0+i*dx);
	  printf("time %.3e PulsarBfield =%.3e\n",pulsarBfield[0][i]/yeartosec,pulsarBfield[1][i]); 
	   
	}
	
	
 
	
}




double getsize(double time){
  double size;
  /*we assume here a free expansion phase (see Van der Swaluw  2001)*/
  C1=pow(6.0/5.0+289/240,-0.2)*pow(source->continuous->L0/source->continuous->pulsartype->E_SN,0.2);
  Vej=sqrt(10.0*source->continuous->pulsartype->E_SN/(3.0*source->continuous->pulsartype->Mej));
 
  size=C1*Vej*pow(time,6.0/5.0);
  printf("sizefunc=%.3e\n",size);
  return size;

}

