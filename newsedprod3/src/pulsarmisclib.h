double getmagneticfield(double time);
double rk4(double(*f)(double, double), double dx, double x, double y);
double Bfield(double Wb, double time);
double rate(double t, double y);
double initialcondition(double time);
void getBfieldlookuptable();
double getsize(double time);
