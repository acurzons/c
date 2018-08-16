struct sourcetype *source;


/*proton distribution array*/

/*proton distribution is a pointer to an array of energy
It will have different size according to what option we ve given into the parameter list:
e.g : option evolution=no -> protondistribution=malloc(sizeof(double*)) 
                             protondistribution[0]=malloc(protonsize*sizeof(double))

     option evolution= yes -> protondistribution=malloc(4*sizeof(double*))
     for (i=0;i>4;i++){
     protondistribution[i]=malloc(protonsize*sizeof(double))
*/

double **protondistribution;
int protonarraysize;
/*electron distribution array*/

double **electrondistribution;
int electronarraysize;
/*The different output from pp interaction*/
struct ARRAYSIZE *gammaraypp;
struct ARRAYSIZE *secondaryelecpp;
struct ARRAYSIZE *neutrinospp;
struct ARRAYSIZE *piondata;
/*The different gamma-ray output from secondary interaction*/
struct ARRAYSIZE *secondaryelecgammasync;
struct ARRAYSIZE *secondaryelecgammabremss;

/*output from electron interaction with photonfield*/
struct ARRAYSIZE *ICgamma;
struct ARRAYSIZE *syncgamma;
struct ARRAYSIZE *bremssgamma;
struct ARRAYSIZE *ICRgamma;   //for IR source 
struct ARRAYSIZE *SSCgamma;

struct ARRAYSIZE *totalgamma;


/*for Pulsar Bfield evolution*/
double *pulsarBfield[2];
double C1; //pulsar free expansion constant 
double Vej; //SNR speed of ejecta
