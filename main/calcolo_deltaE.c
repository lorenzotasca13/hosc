#define MAIN_PROGRAM
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "harmonic_oscillator.h"
#include "random.h"

int main () {

    double **c_tilde, *cbar, **cluster, *deltaE, *error_deltaE;
    int t;


    termalisation();
    c_tilde=C_tilde();
    cbar=Cbar(c_tilde);
    cluster=Cluster(c_tilde, cbar);
    deltaE=DeltaE(cbar);
    error_deltaE=Error_deltaE(deltaE, cluster);

    for (t=0; t<N; t++) printf("%f\t+- %f\n", deltaE[t], error_deltaE[t]);
    

    return 0;
}