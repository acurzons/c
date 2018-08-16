//
//!!**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//!! * bremss_spec.f *                               galprop package * 4/14/2000 
//!!**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//#include "inter.c"
#include "bremss.h"
//#include "simpson.c"






double BremssXS(double Egam, double E0, int IZ1, int Ne1) {
  
  //***********************************************************************
  //          *** I.Moskalenko (MPE,Garching) *** version of 14.05.1998 ***
  // PURPOSE:
  // Calculation of the spectrum of e-bremsstrahlung in the fully ionized,
  // hydrogen-like, or helium-like gas as well as in neutral H, He, and 
  // heavier in Born approximation; valid for E0,Egam >~ 0.01 MeV and
  // beta0, beta1 >> 2*pi*Af*Z, corresponds to gamma=1.001
  // REFERENCES:
  // [1] Blumenthal & Gould 1970,Rev.Mod.Phys.42,237 (BG70); 
  // [2] Gould 1969, Phys.Rev.185,72 (G69);
  // [3] Koch & Motz 1959, Rev.Mod.Phys.31,920 (KM59).
  // [4] Strong A.W., Moskalenko I.V., Reimer O. 2000, ApJ 537, 763
  // INPUT/OUTPUT parameters:
  // Egam, E0 [TeV] are the energy of a gamma and total electron energy,
  // respectively;
  // IZ1 is the nucleus charge; 
  // Ne1 = 0 for unshielded charge, = 1 for 1-e atoms, = 2 for 2-e atoms.
  // dSdK [barn/MeV] is the differential cross section (output).
  //      [now returns  cm^2/TeV   GPR 31/01/2007]
  // Note:
  // #1 a one-parameter Hylleraas function is used for helium-like ions,
  //    while Hartree-Fock functions are used for heutral He atoms;
  // #2 a contribution of free electrons could be estimated by calling
  //    the subroutine with Z=1, Ne=0 and multiplying dSdk by ionization level.
  //***********************************************************************

  Egam *= 1.0e+6; // TeV to MeV GPR
  E0 *= 1.0e+6;   // TeV to MeV GPR

  int i; 
  double eps0,eps1,A,AL,T1,T2,T10,T11,T12,b,OM,phiu,phi1,phi2,p1,p2,Q,dQ;
  double DD1[11];
  double DD[11], PH1[11], PH2[11];  // DD = delta/(2.d0*Af)

  double dSdK = 0.0;
  double mc2 = 0.511;   // MeV, electron rest mass
  double E1 = E0-Egam;  // total energy of an electron after scattering
      
  if (E0 <= 1.02*mc2 || E1 <= mc2) return 1.0e-18*dSdK;
  
  double IZ = IZ1;
  double Ne = Ne1;
  double Pi = 3.141592653590;
  double Af = 7.2974e-3;      // fine structure constant
  double R02 = pow(2.81794e-1,2); // classical electron radius square (barn)
  double EB1 = 0.07; // MeV,kinetic energy; the upper bound for NonRelat.approx.
  double EB2 = 2.0;  // MeV, -"-; the lower bound for HE approx. (nonscreened)
  double delta = Egam*mc2/(2.0*E0*E1);
  double Ef = 0.0;
  double FS = 0.0;
  double FS90 = 0.90;
  double T0 = E0-mc2;   // Ekin
  double P0 = mc2*sqrt(pow((E0/mc2),2)-1.0); // initial electron momentum
  double P1 = mc2*sqrt(pow((E1/mc2),2)-1.0); // final -"-
  double beta0 = P0/E0;
  double beta1 = P1/E1;

  //### Fano-Sauter high-frequency limit (Egam ->E0-mc2, KM59,p.933) ###
  if (Egam/T0 >= FS90) 
    FS =4.*Pi*pow(IZ,3)*R02*Af*Af/T0*beta0*E0/mc2/pow((E0/mc2-1.0),2) 
      *(4./3. +E0/mc2*(E0/mc2-2.0)/(E0/mc2+1.0) *(1.0-log(1.0+2.0*beta0/(1.0-beta0))*pow(mc2/E0,2)/2.0/beta0));
  
  //## Elwert factor (Ef) is valid if   Z*Af*(1./beta1-1./beta0)<<1 ###
  if (T0 < EB2 && (IZ*Af*(1./beta1-1./beta0)) <= 0.5)
    Ef = beta0/beta1 *(1.0-exp(-2.0*Pi*Af*IZ/beta0)) / (1.0-exp(-2.0*Pi*Af*IZ/beta1));
  
  //### NONRELATIVISTIC ENERGIES ( T0 < 0.07 MeV) ###
  //# Born approx.(KM59,p.925,3BNa) with Elwert factor (p.931)
  //# Born approximation: 2*Pi*Z*Af/beta{0,1}<<1; 
  //# nonscreened: 2*delta/[Af*Z^(1/3)] >> 1;
  if(T0 < EB1) {
    A = IZ*IZ;                    // no electron contribution
    dSdK = Ef*R02*Af*A*16./3. /Egam *pow((mc2/P0),2)*log(1.0+2.0*P1/(P0-P1));
    if (dSdK < FS) dSdK = (FS+dSdK)/2.0;   //# 1/2 Fano-Sauter limit
    return 1.0e-18*dSdK;
  }
  
  //### LOW ENEGRY (T0 < 2. MeV) OR delta/(2.*Af*Z) > 4. ###
  //# Born approximation (KM59,p.925,3BN): 2*Pi*Z*Af/beta{0,1}<<1;
  //# nonscreened: 2*delta/Af >> 1
  if (T0 < EB2 || (delta/(2.*Af*IZ)) >= 4.) {
    eps0 = log(1.0+2.0*P0/(E0-P0));  //=dlog((E0+P0)/(E0-P0))
    eps1 = log(1.0+2.0*P1/(E1-P1));  //=dlog((E1+P1)/(E1-P1))
    AL = 2.0*log((E0*E1+P0*P1-mc2*mc2)/(Egam*mc2));
    T1 = 2.0*E0*E1*(pow((P0/P1),2)+1.0)/(P0*P0);
    T2 = mc2*mc2*(eps0*E1/pow(P0,3)+eps1*E0/pow(P1,3)-eps0*eps1/P0/P1);
    T10 = 8.0/3.*(E0*E1)/(P0*P1);
    T11 = Egam*Egam*(pow((E0/P0*E1/P1),2)+1.)/(P0*P1);
    T12 = mc2*mc2 *Egam/(2.0*P0*P1) *(eps0*(E0/P0*E1/P0+1.0)/P0 -eps1*(E0/P1*E1/P1+1.0)/P1+2.0*Egam*E0/P0*E1/P1/(P0*P1));
    A = IZ*IZ + Ne;            // correction for atomic electrons
    if (T0 < EB2) // screening
      A =( IZ*IZ +Ne*(1.-exp(-(T0-EB1)/9./EB1)) ) *(1.-0.3*exp(-Egam/0.033))*Ef;  // left & right connection
    dSdK = R02*Af*A/Egam*P1/P0*(4.0/3.-T1+T2+AL*(T10+T11+T12));
    if(dSdK < FS) dSdK = (FS+dSdK)/2.;   //# 1/2 Fano-Sauter limit
    return 1.0e-18*dSdK;                   
  }
  
  //### HIGH  ENERGY (T0 > 2. MeV) AND delta/(2.*Af*Z) < 4. ###
  //# Schiff approx. for neutral atoms heavier then He (KM59,p.925,3BSe) 
  if(IZ > 2 && IZ == Ne) {
    b = pow(IZ,(1.0/3.0))/111.0/delta;
    OM = 1.0/pow(delta,2)/(1.0+b*b);
    dSdK = 2.0*IZ*IZ*R02*Af/Egam * 
      ((1.0+pow((E1/E0),2)-2.0/3.*E1/E0)*(log(OM)+1.0-2.0/b*atan(b))+
       E1/E0*(2.0/b/b*log(1.0+b*b)+4./3.*(2.0-b*b)/pow(b,3)*atan(b)-8.0/3./b/b+2.0/9.));
    return 1.0e-18*dSdK;
  }
  
  //# arbitrary screening (G69, BG70; KM59,p.925,3BSb)
  A = IZ*IZ + Ne;            //correction for atomic electrons
  phiu = -log(delta)-1.0/2.0;
  //# Hartree-Fock approximation for neutral He atoms
  if (IZ == 2 && IZ == Ne) {
    DD1[1] = log(1.-3);
    for (i=1;i<=10;i++)
      DD1[i] = log(DD[i]);
    
    if (delta/(2.0*Af) >= DD[1]) {
      if(delta/(2.0*Af) <= DD[10]) {
	DINTER(DD1,PH1,11,log(delta/(2.0*Af)),&phi1);
	DINTER(DD1,PH2,11,log(delta/(2.0*Af)),&phi2);
      }
      else {
	phi1 = 4.*A *phiu;   // asymptotics
	phi2 = 4.*A *phiu;
      }
    }
    else {
      phi1 = PH1[1];
      phi2 = PH2[1];
    }
    phi1 = phi1/4.;
    phi2 = phi2/4.;
  }
  else {
    if (Ne == 0) {  //# UNSHIELDED charge
      phi1 = A*phiu;
      phi2 = A*phiu;
    }
    else {
      //c# H-like atoms & Hylleraas-1 approximation for He-like atoms
      if(delta/(2.0*Af) <= 1.0-4) delta = 1.-4*(2.0*Af);
      //SIM1(1.0,delta,1.-3,1.-5,1.-4,FPHI1,*p1);
      //SIM1(1.0,delta,1.-3,1.-5,1.-4,FPHI2,*p2);
      //integrate FPHI1(Q) & FPHI2(Q) from 1.0 to delta in dQ=1.0e-3 steps 
      p1=0.0;
      p2=0.0;
      dQ=1.0e-3;
      for (Q=1.0;Q<=delta;Q+=dQ) {
	p1 = FPHI1(Q, delta, IZ, Ne)*dQ;
	p2 = FPHI2(Q, delta, IZ, Ne)*dQ;
      }

      phi1 = 2.*IZ*(-p1       +1.+(Ne-1.)/IZ)  +pow((IZ-Ne),2)*phiu;
      phi2 = 2.*IZ*(-p2+5./6.*(1.+(Ne-1.)/IZ) )+pow((IZ-Ne),2)*phiu;
    }
  }
  if (phi1 > A*phiu) phi1 = A *phiu;  // asymptotics
  if(phi2 > A*phiu) phi2 = A *phiu;
  
  dSdK = 4.*R02*Af/Egam *((1.0+pow((E1/E0),2))*phi1 -2./3.*E1/E0*phi2);
  if (dSdK < FS) dSdK = (FS+dSdK)/2.;      //`# 1/2 Fano-Sauter limit
  return 1.0e-18*dSdK;
}




double FPHI1(double Q, double delta, int IZ, int Ne) {
  //***********************************************************************
  // used for calculation of PHI1 function (G69,p.74,76)
  //          *** I.Moskalenko (MPE,Garching) *** version of 23.09.1997 ***
  //***********************************************************************
  
  double res = 0.0, AZ, FZ;
  
  double Af = 7.2974e-3;          // fine structure constant
  if (Ne == 1) {      // one-electron atoms
    AZ = 1.0/(2.0*Af*IZ);
    FZ = 1.0/pow((1.0+pow((AZ*Q),2)),2); // form-factor for one-electron atom
    res = pow((Q-delta),2)/pow(Q,3) *(1.0-FZ);
    return res;
  }                   // two-electron atoms (Hylleraas-1 function)
  AZ = 1.0/(2.0*Af*(IZ-5.0/16.0));
  FZ = 1.0/pow((1.0+pow((AZ*Q),2)),2); // form-factor for two-electron atom
  res = pow((Q-delta),2)/pow(Q,3) *(2.*(1.0-FZ) -(1.0-FZ*FZ)/IZ);
  return res;
}

double FPHI2(double Q, double delta, int IZ, int Ne) { 
  //***********************************************************************
  // used for calculation of PHI2 function (G69,p.74,76)
  //          *** I.Moskalenko (MPE,Garching) *** version of 23.09.1997 ***
  //***********************************************************************
    
  double res = 0.0, AZ, FZ;

  double Af = 7.2974e-3;          // fine structure constant
  if (Ne == 1) {      // one-electron atoms
    AZ = 1.0/(2.0*Af*IZ);
    FZ = 1.0/pow((1.0+pow((AZ*Q),2)),2); // atomic form-factor
    res=(pow(Q,3)+(3.0-6.0*log(Q/delta))*delta*delta*Q-4.0*pow(delta,3))/pow(Q,4)*(1.0-FZ);
    return res;
  }                        // two-electron atoms (Hylleraas-1 function)
  AZ = 1.0/(2.0*Af*(IZ-5.0/16.0));
  FZ = 1.0/pow((1.0+pow((AZ*Q),2)),2); // form-factor for two-electron atom
  res=(pow(Q,3)+(3.0-6.0*log(Q/delta))*delta*delta*Q-4.0*pow(delta,3))/pow(Q,4) *(2.*(1.0-FZ) -(1.0-FZ*FZ)/IZ);
  return res;
}
