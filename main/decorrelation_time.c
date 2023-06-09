#define MAIN_PROGRAM
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "harmonic_oscillator.h"
#include "random.h"

int main () {

    int t_M;
    double** gamma;

    FILE *fd;
    fd=fopen("gamma.txt", "w");
   
    termalisation();

    gamma=gamma_decorr();
    for (t_M=0; t_M<t_M_max; t_M++) {
        fprintf(fd, "%d\t%f\n", t_M, gamma[3][t_M]/gamma[3][0]);
    }

    fclose(fd);

    return 0;
}