/*library for the protondistribution functions*/
double *getgamma1andGMINcontinuousproton();
double *protonGMINGMAX();
void constantprotondistribution();
void diracprotondistribution();
double normalisationprotoncontinuous(double time);
double  pulsarenergy (double time);
double normalisationprotonimpulsive();
double normalisationprotondirac();
void evolutionprotoncontinuous();
double getinitialcontinuousprotondistribution(double gamma,double tau);
double getescapeprotoncontinuous(double gamma,double tau);
void evolutionprotonimpulsive();
double getinitialimpulsiveprotondistribution(double gamma);
