#define MAIN_PROGRAM
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "harmonic_oscillator.h"
#include "random.h"

int main () {

    double** ctilde;
    int t,k;

    termalisation();
    ctilde=C_tilde();
    for (t=0; t<N; t++) {
		for (k=0; k<N_BIN; k++) {
        printf("%f\n", ctilde[t][k]);
        }
	}
    

    
    return 0;
}