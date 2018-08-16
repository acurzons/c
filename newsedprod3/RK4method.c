#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define eta 0.03
#define L0 3.1E39 //erg/s
#define E0 1E51 //ergs
#define solarmass  1.9891E33 //in g
#define Mej (9.5*solarmass) // in solarmass
#define yeartosec 3.156e+7 //sec
#define t0 (730*yeartosec) //yr
#define brakingindex 2.509 

double C1;
double Vej;

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
  return eta*L0*pow(1.0+t/t0,-(brakingindex+1.0)/(brakingindex-1.0))-6.0/5.0*y/t;
  //return eta*L0-6.0/5.0*y/t;

}

double initialcondition(double time){
 
  return 5.0/11.0*eta*L0*time; 



}
int main(void)
{
	double *y, x, y2;
	double x0 = 5*yeartosec, x1 = 10000*yeartosec, dx =0.5*yeartosec;
	int i, n = 1 + (x1 - x0)/dx;
	FILE *ofile;
	double B,B2;
	double source;
	C1=pow(6.0/5.0+289/240,-0.2)*pow(L0/E0,0.2);
	Vej=sqrt(10.0*E0/(3.0*Mej));
	printf("C1=%.3e  Vej=%.3e\n",C1,Vej);
	sleep(10);
	if((ofile=fopen("testRK4.txt","w+"))==NULL){
	  printf("unable to open file \n");
	  exit(1);
	}
	y = malloc(sizeof(double) * n);
 
	for (y[0] = initialcondition(x0), i = 1; i < n; i++){
	  y[i] = rk4(rate, dx, x0 + dx * (i - 1), y[i-1]);
	  B=Bfield(y[i],x0+i*dx);
	 
	  B2=sqrt(30.0/11.0*eta*L0/(pow(C1*Vej,3.0)))*pow(x0+i*dx,-1.3);
	  source=eta*L0*pow(1.0+(x0+i*dx)/t0,-(brakingindex+1.0)/(brakingindex-1.0));
	  printf("y[%d]=%.3e %.3e %.3e %.3e %.3e %.3e\n",i,y[i],B,source,(x0+i*dx)/t0,t0,x0+i*dx);
	  // sleep(1);
	  fprintf(ofile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",(x0+i*dx)/yeartosec,y[i],B,B2,source,eta*L0);
	}
	fclose(ofile);
 
	printf("x\ty\trel. err.\n--------x0=%.3e t0=%.3e\n",x0,t0);
	
 
	return 0;
}
