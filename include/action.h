#ifndef HOSC_H
#define HOSC_H
#endif

#ifndef ACTION_C
extern double action();
extern double delta_action(int , double);
extern double sweep (int);
extern double* corr(int t_fisico);
extern double** gamma_decorr();
extern void termalisation();
extern double c(int t_fisico);
extern double** C_tilde();
extern double* Cbar(double**);
extern double ** Cluster(double**, double*);
extern double arcosh(double x);
extern double* DeltaE(double*);
extern double* Error_deltaE(double*, double**);

#endif
