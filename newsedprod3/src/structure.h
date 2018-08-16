/* this file host the struct variable declaration */

struct powerlaw{
  double inp1;
  double GMIN;
  double GMAX;
  double EMIN;
  double EMAX;
};

struct brokenpowerlaw{
  double inp1;
  double inp2;
  double GMIN;
  double GMAX;
  double Ebreak;
  double EMIN;
  double EMAX;
};

struct powerlawcutoff{
  double inp1;
  double expinp;
  double GMIN;
  double GMAX;
  double Ecutoff;
  double EMIN;
  double EMAX;
};

struct brokenpowerlawcutoff{
  double inp1;
  double inp2;
  double expinp;
  double GMIN;
  double GMAX;
  double Ebreak;
  double Ecutoff;
  double EMIN;
  double EMAX;
};


struct impulsive{
  double E0;
  char spectratypeelectron[100];
  char spectratypeproton[100];
  char evolution[10];
  char escape[10];
  double Efractionproton;
  double distance;
  double size;
  double density;
  double Bfield;
  double time;
  double TempIR;
  double EnergyIR;
  struct powerlaw *powerlawelectron;
  struct brokenpowerlaw *brokenpowerlawelectron;
  struct powerlawcutoff *powerlawcutoffelectron;
  struct brokenpowerlawcutoff *brokenpowerlawcutoffelectron;
  struct powerlaw *powerlawproton;
  struct brokenpowerlaw *brokenpowerlawproton;
  struct powerlawcutoff *powerlawcutoffproton;
  struct brokenpowerlawcutoff *brokenpowerlawcutoffproton;
};

struct pulsartype{
  double brakingindex;
  double period; //in s
  double perioddot; //in s s-1
  double fractionmag; 
  double Mej; //the input is in solar mass but will be converted in g
  double E_SN;
  double t0; //critical value
};

struct continuous{
  double L0;
  char spectratypeelectron[100];
  char spectratypeproton[100];
  char evolution[10];
  char escape[10];
  double Efractionproton;
  double distance;
  double size;
  double density;
  double Bfield;
  double time;
  double TempIR;
  double EnergyIR;
  char pulsar[10];
  struct powerlaw *powerlawelectron;
  struct brokenpowerlaw *brokenpowerlawelectron;
  struct powerlawcutoff *powerlawcutoffelectron;
  struct brokenpowerlawcutoff *brokenpowerlawcutoffelectron;
  struct powerlaw *powerlawproton;
  struct brokenpowerlaw *brokenpowerlawproton;
  struct powerlawcutoff *powerlawcutoffproton;
  struct brokenpowerlawcutoff *brokenpowerlawcutoffproton;
  struct pulsartype *pulsartype;
};

struct dirac{
  double E0;
  double Edirac;
  double Bfield;
  double size;
  double density;
  double distance;
  double TempIR;
  double EnergyIR;
};

struct sourcetype {
  char sourcetype[50];
  struct dirac *dirac;
  struct impulsive *impulsive;
  struct continuous *continuous;
};

struct ARRAYSIZE {
  double *array[2];
  int size;
};
