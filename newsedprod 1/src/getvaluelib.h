void getlistparameters();
struct powerlaw *getspectrumpowerlaw(FILE *infile,char filename[200]);
struct brokenpowerlaw *getspectrumbrokenpowerlaw(FILE *infile,char filename[200]);
struct powerlawcutoff *getspectrumpowerlawcutoff(FILE *infile,char filename[200]);
struct brokenpowerlawcutoff *getspectrumbrokenpowerlawcutoff(FILE *infile,char filename[200]);
struct pulsartype *getpulsarparameter(FILE *infile, char filename[200]);
