#define MAIN_PROGRAM
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "harmonic_oscillator.h"
#include "random.h"

int main () {
    int i;
    int heat=0; /*0 a freddo, 1 a caldo*/
    double accettanza_media=0;
    int n_sweep=200; /*numero di sweep effettuati*/
    FILE *fd;
    fd=fopen("term.txt", "w");

    if (heat==0) { /*inizializziamo il vettore*/
        for (i=0; i<N; i++) {
            xx[i]=0;
        }
    } 
    else {
        rlxd_init(1,3122000); 
	    ranlxd(xx,N);
    }
    
    for (i=0; i<n_sweep; i++) {
        fprintf(fd, "%d\t%f\n", i, action());
        accettanza_media+=sweep(3122000+i);
    
    }
    fclose(fd);
    accettanza_media/=(double)n_sweep;
    printf("l'accettanza media è: %f\n",accettanza_media);
    return 0;
}