#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structure.h"
#include "constant.h"
#include <string.h>
#include "globalparameters.h"
#include "protondistributionlib.h"
#include "protonSEDlib.h"

void ppinteractiongamma(){
  /*This functions aims at getting the main secondary electrons+neutrinos+gamma-ray from the given proton energy distribution, we remember that the dynamic array start with GMIN or gamma2 (see Manolakou et al 2007 for gamma2), thus we must expect */

  /*We are reminded that the main  energy unit here is TeV*/

  /*First the gammaloop*/
  int i,j; //loop integer variable
  double dEgamma;
  double dEp,Ep;
  double Tp,dTp;
  double x,dx;
  double *GMINMAX;
  double EGMIN;
  int sizegamma;
  /* for the secondary electron*/
  
  int sizesecelectron;
  
  double dEsec;

  /*for the neutrino*/
  int sizeneutrino;
 double dEneut;

  double density;
  if(source->continuous!=NULL) density=source->continuous->density;
  else if(source->impulsive!=NULL) density=source->impulsive->density;
  GMINMAX=protonGMINGMAX();
  EGMIN=GMINMAX[0]*protonmassTeV;
 
  /*initialisation of the gammaraypp*/
  sizegamma=(int)(20+8)*100/STEP;
  gammaraypp=malloc(sizeof(struct ARRAYSIZE));
  gammaraypp->size=sizegamma;
  for(i=0;i<2;i++){
    gammaraypp->array[i]=malloc(sizegamma*sizeof(double));
  }
  for(i=0;i<2;i++){
    for(j=0;j<sizegamma;j++){
      gammaraypp->array[i][j]=0;
    }
  }
  /*First the gamma-ray from proton-proton innteraction*/

  
  for(i=0;i<=sizegamma;i++){
    gammaraypp->array[0][i]=pow(10,(double)(i*STEP-2000.0)/100.0);

    dEgamma=gammaraypp->array[0][i]*(pow(10,(double)STEP/100.0)-1);
    
    
    
    for(j=0;j<protonarraysize;j++){
 
      Ep=protondistribution[0][j]*protonmassTeV;
      dEp=Ep*(pow(10,(double)STEP/100.0)-1); //separation between the 2 energy in array
      Tp=Ep-protonmassTeV;
      dTp=Tp*(pow(10,(double)STEP/100.0)-1);
      
      if ((Tp>0.2797E-3)&&(gammaraypp->array[0][i]<Ep)) gammaraypp->array[1][i] +=lightspeedcm*density*protondistribution[1][j]*getAmax(Ep)*newFgamma(gammaraypp->array[0][i],Ep)*dEp;
      

       

    }

  }

  
  sizesecelectron=(int)(20+8)*100/STEP+1;

  secondaryelecpp=malloc(sizeof(struct ARRAYSIZE));
  secondaryelecpp->size=sizesecelectron;
  for(i=0;i<2;i++){
    secondaryelecpp->array[i]=malloc(sizesecelectron*sizeof(double));
  }
  for(i=0;i<2;i++){
    for(j=0;j<sizegamma;j++){
      secondaryelecpp->array[i][j]=0;
    }
  }

  /*The secondary electron coming from proton-proton interaction*/

  
  for(i=0;i<=sizesecelectron;i++){
    secondaryelecpp->array[0][i]=pow(10,(double)(i*STEP-2000)/100.0);

    dEsec=secondaryelecpp->array[0][i]*(pow(10,(double)STEP/100.0)-1);
    
    
    
    
    for(j=0;j<protonarraysize;j++){
      Ep=protondistribution[0][j]*protonmassTeV;
      dEp=Ep*(pow(10,(double)STEP/100.0)-1); //separation between the 2 energy in array

   
      if((Ep>=secondaryelecpp->array[0][i])&&(secondaryelecpp->array[0][i]>electronmassTeV)){
	secondaryelecpp->array[1][i] +=lightspeedcm*density*protondistribution[1][j]*PPCrossSection(Ep)*Felectron(secondaryelecpp->array[0][i]/Ep,Ep)*dEp/Ep;

      }

    }

  }

  /*This part with the neutrino flux will need to be tested thoroughly, I need to look more thoroughly at the method needed to compute the different functions */
  sizeneutrino=28*100.0/STEP+1;
  neutrinospp=malloc(sizeof(struct ARRAYSIZE));
  neutrinospp->size=sizeneutrino;
  for(i=0;i<2;i++){
    neutrinospp->array[i]=malloc(sizeneutrino*sizeof(double));
  }
  for(i=0;i<2;i++){
    for(j=0;j<sizegamma;j++){
      neutrinospp->array[i][j]=0;
    }
  }





for(i=0;i<=sizeneutrino;i++){
    neutrinospp->array[0][i]=pow(10,(double)(i*STEP-8000.0)/100.0);

    dEneut=gammaraypp->array[0][i]*(pow(10,(double)STEP/100.0)-1);
    
    
    
    
    for(j=0;j<protonarraysize;j++){
      Ep=protondistribution[0][j]*protonmassTeV;
      dEp=Ep*(pow(10,(double)STEP/100.0)-1); //separation between the 2 energy in array
      if((Ep>=neutrinospp->array[0][i])&&(neutrinospp->array[0][i]>neutrinomassTeV)){
	neutrinospp->array[1][i] +=lightspeedcm*density*PPCrossSection(Ep)*protondistribution[1][j]*(Fgamma(neutrinospp->array[0][i]/Ep,Ep)+Fneutrino(neutrinospp->array[0][i]/Ep,Ep))*dEp/Ep;

      }

    }

 }

}


double Fgamma(double x, double Ep) 
{ // Fgamma function Eq 58 of Kelner etal 2006

  double L, Bg, betag, kg, Fg, Fg1, Fg2, Fg3, Fg4, xbg;

  L = log(Ep);
  Bg = 1.30+0.14*L+0.011*L*L;
  betag = 1.0 / (1.79+0.11*L+0.008*L*L);
  kg = 1.0 / (0.801+0.049*L+0.014*L*L);

  xbg = pow(x, betag);

  Fg1 = (1.0-xbg)/(1.0+kg*xbg*(1.0-xbg));
  Fg2 = 1.0/log(x);
  Fg3 = 4.0*betag*xbg/(1.0-xbg); 
  Fg4 = 4.0*kg*betag*xbg*(1.0-2.0*xbg)/(1.0+kg*xbg*(1.0-xbg));
  Fg  = Bg*log(x)/x*pow(Fg1,4.0)*(Fg2-Fg3-Fg4);

  return Fg;
}



double PPCrossSection(double Ep)
{ // inclusive pp-->2gamma cross section Eq 79 Kelner etal 2006

  double xs = 0.0;
 
  double Eth = 1.0e-12 * (protonmassTeV*1E12 + 2.0*pionmassTeV*1E12 + pow(pionmassTeV*1E12,2)/(2.0*protonmassTeV*1E12));
  double Tp=Ep-protonmassTeV;
  double Tth=2.797E-4; //TeV
  double L = log(Tp/Tth);
  //double L=log(Ep/1.0);
  // xs = (34.3+1.88*L+0.25*L*L) * pow((1.0-pow((Eth/Ep),4.0)),2.0); // units mb
  xs = (30.7-0.96*L+0.18*L*L) * pow((1.0-pow((Tth/Tp),1.9)),3.0); // units mb
  xs *= 1.0e-27; // units cm^2

  if (Ep <= Eth){
     xs = 0.0;
     // printf("xs=%lf\n",xs);
  }
  //printf("  PPCrossSection %.6lf %lf %lf %e\n", Ep, L, Eth, xs);

  return xs;
}



double Felectron(double x, double Ep) 
{ // Felectron function Eq 62 of Kelner etal 2006
  // This is also valid for muon-decay neutrinos

  double L = log(Ep);
  double Be = 1.0/(69.5+2.65*L+0.3*L*L);
  double betae = 1.0 / pow(0.201+0.062*L+0.00042*L*L,0.25);
  double ke = (0.279+0.141*L+0.0172*L*L)/(0.3+pow(2.3+L,2.0));

  double Fe, xbe;

  xbe = pow(x, betae);
 
  Fe = Be * (pow(-1.0*log(x),5.0)) * pow(1.0+ke*pow(log(x),2.0),3.0)/(x*(1.0+0.3/xbe));
  if(Fe!=Fe) Fe=0.0;
  return Fe;
}


double Fneutrino(double x, double Ep) 
{ // Fneutrino function Eq 66 of Kelner etal 2006
  // for pion --> muon + neutrino

  double L = log(Ep);
  double Bn = 1.75+0.204*L+0.010*L*L;
  double betan = 1.0 / (1.67+0.111*L+0.0038*L*L);
  double kn = 1.07-0.086*L+0.002*L*L;
  double y = x/0.427;

  double Fn, xbn;

  if (x > 0.427) return 0.0;

  xbn = pow(y, betan);

  Fn = Bn * log(y)/y * pow((1.0-xbn)/(1.0+kn*xbn*(1.0-xbn)),4.0) * 
    (1.0/log(y) - (4.0*betan*xbn)/(1.0-xbn) - (4.0*kn*xbn*(1.0-2.0*xbn))/(1.0+kn+xbn*(1.0-xbn)));
  
  return Fn;
}




void ppinteractiondirac(){
  /*This functions aims solely give the neutral  pion gamma-ray production from a monoenergetic set of protons*/

  /*We are reminded that the main  energy unit here is TeV*/

  /*First the gammaloop*/
  int i,j; //loop integer variable
  double dEgamma;
  double dEp,Ep;
  double *GMINMAX;
  double EGMIN;
  int sizegamma;
  /* for the secondary electron*/
  
  int sizesecelectron;
  
  double dEsec;

  /*for the neutrino*/
  int sizeneutrino;
 double dEneut;

  double density;
  density=source->dirac->density;
  
  GMINMAX=protonGMINGMAX();
  EGMIN=GMINMAX[0]*protonmassTeV;
 
  /*initialisation of the gammaraypp*/
  sizegamma=(int)(20+8)*100/STEP;
  gammaraypp=malloc(sizeof(struct ARRAYSIZE));
  gammaraypp->size=sizegamma;
  for(i=0;i<2;i++){
    gammaraypp->array[i]=malloc(sizegamma*sizeof(double));
  }
  for(i=0;i<2;i++){
    for(j=0;j<sizegamma;j++){
      gammaraypp->array[i][j]=0;
    }
  }
  /*First the gamma-ray from proton-proton innteraction*/

  
  for(i=0;i<=sizegamma;i++){
    gammaraypp->array[0][i]=pow(10,(double)(i*STEP-2000.0)/100.0);

    dEgamma=gammaraypp->array[0][i]*(pow(10,(double)STEP/100.0)-1);
    
    
    
    
    
      Ep=protondistribution[0][0]*protonmassTeV;
      dEp=Ep*(pow(10,(double)STEP/100.0)-1); //separation between the 2 energy in array
      if(Ep>=gammaraypp->array[0][i]){
	gammaraypp->array[1][i] +=4*PI*lightspeedcm*density*newFgamma(gammaraypp->array[0][i],Ep)*getAmax(Ep)*protondistribution[1][0];

      }

    

  }
printf("who let the dogs out\n");
}


double interpolationprotonflux(double protonenergy){
  int i;
  double lowerbound, lowerboundflux;
  double upperbound, upperboundflux;
  double interpolatedflux;
  /* WE HAVE TO MAKE SURE THE PROGRAM NEVER GO BEYOND THE MAXIMUM ARRAY*/
  if ((protonenergy<protondistribution[0][0]*protonmassTeV)&&(protonenergy>protondistribution[0][protonarraysize-1]*protonmassTeV)){
    interpolatedflux=0.0;
    return interpolatedflux;
  }

  while (protondistribution[0][i]*protonmassTeV<protonenergy){
    i++;
  }
  /*Now we have the upper 'i' and lower bound 'i-1'*/
  lowerbound=protondistribution[0][i-1]*protonmassTeV;
  upperbound=protondistribution[0][i]*protonmassTeV;
  lowerboundflux=protondistribution[1][i-1];
  upperboundflux=protondistribution[1][i];
  
  interpolatedflux=(log10(upperboundflux)-log10(lowerboundflux))/(log10(upperbound)-log10(lowerbound))*(log10(protonenergy)-log10(lowerbound))+log10(lowerboundflux);
  interpolatedflux=pow(10,interpolatedflux);




  return interpolatedflux;
}
  
  
double pionXS(double Ep){
  /*return the cross section  for the production of pion*/
  const double pionmassGeV=pionmassTeV*1E3;
  const double protonmassGeV=protonmassTeV*1E3;
  const double Mres=1.1883; //GeV
  const double Gres=0.2264; //GeV
  const double Tthres=0.2797; // KINETIC energy threshold in GeV
  const double sigma0=7.66E-30; //cm^2
  double pioncrosssection;
  double Tp=Ep*1E3-protonmassGeV; //kinetic energy in GeV
  double sigmainel=(30.7-0.96*log(Tp/Tthres)+0.18*pow(log(Tp/Tthres),2.0))*pow(1.0-pow(Tthres/Tp,1.9),3.0)*1E-27; //inelastic XS in cm^2
  
  double s=2.0*protonmassGeV*(Tp+2.0*protonmassGeV);
  double eta=sqrt(pow(s-pow(pionmassGeV,2.0)-4.0*pow(protonmassGeV,2.0),2.0)-16.0*pow(pionmassGeV,2.0)*pow(protonmassGeV,2.0))/(2.0*pionmassGeV*sqrt(s));
  double gamma=sqrt(pow(Mres,2.0)*(pow(Mres,2.0)+pow(Gres,2.0)));
  double K=sqrt(8.0)*Mres*Gres*gamma/(PI*sqrt(pow(Mres,2.0)+gamma));
  double fbw=protonmassGeV*K/(pow(pow(sqrt(s)-protonmassGeV,2.0)-pow(Mres,2.0),2.0)+pow(Mres*Gres,2.0));
  
  double sigma1p=sigma0*pow(eta,1.95)*(1.0+eta+pow(eta,5.0))*pow(fbw,1.86);
  double sigma2p;
  if (Tp>0.56) sigma2p=5.7E-27/(1.0+exp(-9.3*(Tp-1.4)));
  else  sigma2p=0.0;
  if(Tp<2.0) pioncrosssection=sigma1p+sigma2p;
  else pioncrosssection=sigmainel*averagepionmultiplicity(Tp);

  return pioncrosssection;
    

}
 

double averagepionmultiplicity(double Tp){
  /* we will use the GEANT MODEL from the latest parametrisation method*/
  const double Tthres=0.2797; // KINETIC energy threshold in GeV
  const double protonmassGeV=protonmassTeV*1E3;
  double Qp=(Tp-Tthres)/(protonmassGeV);
  double pionmultiplicity;

  if((Tp>=1.0)&&(Tp<5.0)){
    pionmultiplicity=-6E-3+0.237*Qp-0.023*pow(Qp,2.0);
  }

  else if(Tp>=5.0){
    double a1=0.728,a2=0.596,a3=0.491,a4=0.2503,a5=0.117;
    double ksip=(Tp-3.0)/(protonmassGeV);
    if ((QGSJET==1)&&(Tp>100.0)) {
      a1=0.908;
      a2=0.0009;
      a3=6.089;
      a4=0.176;
      a5=0.448;
    }
    pionmultiplicity=a1*pow(ksip,a4)*(1.0+exp(-a2*pow(ksip,a5)))*(1.0-exp(-a3*pow(ksip,1.0/4.0)));
  }

  return pionmultiplicity;



}   


double getAmax(double Ep){
  double Amax;
  const double pionmassGeV=pionmassTeV*1E3;
  const double protonmassGeV=protonmassTeV*1E3;
  const double Tthres=0.2797;
  double Tp=Ep*1E3-protonmassGeV; //kinetic energy in GeV
  double s=2.0*protonmassGeV*(Tp+2.0*protonmassGeV);

  double EpionCM=(s-4*pow(protonmassGeV,2.0)+pow(pionmassGeV,2.0))/(2*sqrt(s));
  double gammaCM=(Tp+2*protonmassGeV)/(sqrt(s));
  double PpionCM=sqrt(pow(EpionCM,2.0)-pow(pionmassGeV,2.0));;
  double betapionCM=sqrt(1.0-pow(gammaCM,-2.0));
  double Emaxpionlab=gammaCM*(EpionCM+PpionCM*betapionCM);
  double thetap=Tp/protonmassGeV;

  /*Now let s compute the Amax value */

  if(Tp<1.0) Amax=5.9*pionXS(Ep)/(Emaxpionlab);
  else { /* For Giant there are different fit for different energy band*/
    double b1,b2,b3;
    if((Tp>=1.0)&&(Tp<5.0)){
      b1=9.53;
      b2=0.52;
      b3=0.054;
    }
    else{
      b1=9.13;
      b2=0.35;
      b3=9.7E-3;
      if((QGSJET==1)&&(Tp>100)){
	b1=13.16;
	b2=0.4419;
	b3=0.01439;
      }
    }
    Amax=b1*pow(thetap,-b2)*exp(b3*pow(log(thetap),2.0))*pionXS(Ep)/protonmassGeV;

  }
  return Amax*1E3; //in cm^2 TeV^-1
}



double mu(double  Tp){
  double result;
    const double protonmassGeV=protonmassTeV*1E3;
    double q=(Tp-1.0)/protonmassGeV;  
    result=5.0/4.0*pow(q,5.0/4.0)*exp(-5.0/4.0*q);
    return result;
  }

double kappa(double Tp){
  double k;
  const double protonmassGeV=protonmassTeV*1E3;
  double thetap=Tp/protonmassGeV;
  k=3.29-1.0/5.0*pow(thetap,-3.0/2.0);
  return k;
}

double newFgamma(double Egamma,double Ep) {
  double fgamma;
  Egamma *=1.0E3; // in GeV
  const double pionmassGeV=pionmassTeV*1E3;
  const double protonmassGeV=protonmassTeV*1E3;
  const double Tthres=0.2797; // KINETIC energy threshold in GeV
  double Tp=Ep*1E3-protonmassGeV; //kinetic energy in GeV
  double s=2.0*protonmassGeV*(Tp+2.0*protonmassGeV);


  double EpionCM=(s-4*pow(protonmassGeV,2.0)+pow(pionmassGeV,2.0))/(2*sqrt(s));
  double gammaCM=(Tp+2*protonmassGeV)/(sqrt(s));
  double PpionCM=sqrt(pow(EpionCM,2.0)-pow(pionmassGeV,2.0));
  double betapionCM=sqrt(1.0-pow(gammaCM,-2.0));
  double Emaxpionlab=gammaCM*(EpionCM+PpionCM*betapionCM);
  double gammapionlab=Emaxpionlab/pionmassGeV;
  double betapionlab=sqrt(1.0-pow(gammapionlab,-2.0));
  double Egammamax=pionmassGeV/2.0*gammapionlab*(1.0+betapionlab);
  double thetap=Tp/protonmassGeV;
  

  /*Fractions*/
  double Yy=Egamma+pow(pionmassGeV,2.0)/(4.0*Egamma);
  double Yymax=Egammamax+pow(pionmassGeV,2.0)/(4.0*Egammamax);
  double Xy=(Yy-pionmassGeV)/(Yymax-pionmassGeV);

  double lambda,alpha,beta,gamma;
  if (Tp<1.0){
    lambda=1.0;
    alpha=1.0;
    beta=kappa(Tp);
    gamma=0.0;
  }

  else if ((Tp>=1.0)&&(Tp<=4.0)){
    lambda=3.0;
    alpha=1.0;
    beta=mu(Tp)+2.45;
    gamma=mu(Tp)+1.45;
  }

  else if ((Tp>4.0)&&(Tp<=20.0)){
    lambda=3.0;
    alpha=1.0;
    beta=3.0/2.0*mu(Tp)+4.95;
    gamma=mu(Tp)+1.50;
  }
      
  else if((Tp>20.0)&&(Tp<=100.0)){
    lambda=3.0;
    alpha=0.5;
    beta=4.2;
    gamma=1.0;
  }
  else if (Tp>100.0){
    lambda=3.0;
    alpha=0.5;
    beta=4.9;
    gamma=1.0;
    if (QGSJET==1){
      lambda=3.55;
      alpha=0.5;
      beta=4.5;
      gamma=1.0;
    }
  }
  double C=lambda*pionmassGeV/Yymax;
  fgamma=pow(1.0-pow(Xy,alpha),beta)/(pow(1.0+Xy/C,gamma));
    if(fgamma!=fgamma){
      fgamma=0.0;
    }


return fgamma;
}
