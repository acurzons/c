/*library for the electrondistribution functions*/
double *getgamma1andGMINcontinuouselectron();
double *electronGMINGMAX();
void constantelectrondistribution();
void diracelectrondistribution();
double normalisationelectroncontinuous(double time);
double  pulsarenergy (double time);
double normalisationelectronimpulsive();
double normalisationelectrondirac();
void evolutionelectroncontinuous();
double getinitialcontinuouselectrondistribution(double gamma,double tau);
double getescapeelectroncontinuous(double gamma,double tau);
void evolutionelectronimpulsive();
double getinitialimpulsiveelectrondistribution(double gamma);
