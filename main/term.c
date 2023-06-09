#define MAIN_PROGRAM
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "harmonic_oscillator.h"
#include "random.h"

int main () {
    int i;
    int heat=1; /*0 a freddo, 1 a caldo*/
    double accettanza_media=0;
    FILE *fd;
    fd=fopen("termalisation_hot.txt", "w");

    if (heat==0) { /*inizializziamo il vettore*/
        for (i=0; i<N; i++) {
            xx[i]=0;
        }
    } 
    else {
        rlxd_init(1,3122000); 
	    ranlxd(xx,N);
        for (i=0; i<N; i++) xx[i]*=2;
    }
    
    for (i=0; i<N_TERM; i++) {
        fprintf(fd, "%d\t%f\n", i, action());
        accettanza_media+=sweep(3122000+i);
    
    }
    fclose(fd);
    accettanza_media/=(double)N_TERM;
    printf("L'accettanza media è: %f\n",accettanza_media);
    return 0;
}