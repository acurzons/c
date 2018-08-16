#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define yeartosec 3.156e+7 //sec


int main(int argc,char **argv){
  double P0,P;
  double Pdot0,Pdot;
  double brakingindex;
  double time;
  P0=atof(argv[1])*1E-3; //input in ms
  Pdot0=atof(argv[2]); // input in sec
  brakingindex=atof(argv[3]);
  time=atof(argv[4])*yeartosec;
  P=P0*pow(1+((brakingindex-1.0)*Pdot0*time)/P0,1.0/(brakingindex-1.0));
  Pdot=pow(P0/P,brakingindex-2.0)*Pdot0;

  printf("-------------------------\n");
  printf("P(%.1f)=%.3f ms and dotP=%.3e s s-1\n",time/yeartosec,P*1E3,Pdot);
printf("-------------------------\n");


  return 0;
}
