#define ACTION_C
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "random.h"
/*in generale k è il tempo fisico mentre i quello markoviano*/

double action () {
	int i;
	double S=0;
	for (i=0; i<N; i++) {
		S=S+(xx[(i+1)%N]-xx[i])*(xx[(i+1)%N]-xx[i])+W*W*xx[i]*xx[i];
	}
	S=S*M/2;
	return S;

}
	

double delta_action (int k, double deltax) {

	double deltaS=deltax*M/2*(-2*xx[(k+1)%N]+(2*W*W+4)*xx[k]-2*xx[(k-1+N)%N]+2*deltax+W*W*deltax);
	return deltaS;

}

double sweep (int seed) { 
	
	double r[2*N];	
	double deltaS;
	double deltax=0;
	int i;
	double accettanza=0;

	rlxd_init(1,seed);
	ranlxd(r,2*N);
	for (i=0; i<N; i++) {
		deltax=2*Delta*(r[i]-0.5);
		deltaS=delta_action(i, deltax);
		if (deltaS<0 || r[N+i]<=exp(-deltaS)) {
			xx[i]=xx[i]+deltax;
			accettanza+=1.0;
		}
	}
	accettanza/=(double)N;
	return accettanza;
}

/*inizializza il vettore e lo termalizza*/
void termalisation() {

	int i;
    int heat=0; /*0 a freddo, 1 a caldo*/
	int n_sweep=200; /*numero di sweep effettuati*/

	if (heat==0) { /*inizializziamo il vettore*/
        for (i=0; i<N; i++) {
            xx[i]=0;
        }
    } 
    else {
        rlxd_init(1,3122000); 
	    ranlxd(xx,N);
    }

	for (i=0; i<n_sweep; i++) { /*fa un po' di sweep per termalizzare*/
		sweep(3122000+i);
    }


	return;
}

/*corr[i] è la correlazione sull'i-esimo cammino, a t_fisico fissato
double* corr(int t_fisico) {

	double* corr=malloc(N_SWEEP*sizeof(double));
	int k,i;
	
	for (i=0; i<N_SWEEP; i++) {
		corr[i]=0;
		sweep(15092000+i);
		for (k=0; k<N; k++) {
			corr[i]+=xx[k]*xx[(k+t_fisico)%N];
		}
		corr[i]/=(double)N;
	}

	return corr;
}
*/
/*calcola c(t) sul xx attuale*/
double c(int t_fisico) {
	double c=0;
	int k;

	for (k=0; k<N; k++) {
		c+=xx[k]*xx[(k+t_fisico)%N];	
	}
	c/=(double)N;

	return c;
}

/*calcola gamma, che è funzione di t_M e t_fisico*/
/*double** gamma_decorr() {

	int i,t_fisico, t_M;
	double** gamma=malloc(t_fisico_max*sizeof(double*));
	double* correl;
	double y_medio;
	for (i=0; i<t_fisico_max; i++) gamma[i]=malloc(t_M_max*sizeof(double));
	
	
	for (t_fisico=0; t_fisico<t_fisico_max; t_fisico++) {
		for (t_M=0; t_M<t_M_max; t_M++) {
			sweep(t_M+34567);
			
		}
	}




	for (t_fisico=0; t_fisico<t_fisico_max; t_fisico++) {
		correl=corr(t_fisico);
		for (t_M=0; t_M<t_M_max; t_M++) {
			gamma[t_fisico][t_M]=0;
			y_medio=0;
			for (i=0; i<N_SWEEP-15; i++) {
				gamma[t_fisico][t_M]+=c(i)*c(i+t_M);
				y_medio+=correl[i];
			}
			y_medio/=(double)(N_SWEEP-15);
			gamma[t_fisico][t_M]/=(double)(N_SWEEP-15);
			gamma[t_fisico][t_M]-=y_medio*y_medio;
		}
	}
	
	return gamma;
}

restituisce un vettore lungo N_bin con le c_tilde, ovvero le medie sui bin delle c(t)*/



double** C_tilde() { /*C_tilde[t_fisico][t_mark]*/

	int t,k;
	double **c_tilde,*p;
	c_tilde=(double**)malloc(N*sizeof(double*));
	p=(double*)malloc(N_BIN*N*sizeof(double));
	for (t=0; t<N; t++) c_tilde[t]=p+t*N_BIN;
	for (t=0; t<N; t++) {
		for (k=0; k<N_BIN; k++) c_tilde[t][k]=0;
	}
    /*fine dichiarazioni*/

	for (k=0; k<N_SWEEP; k++) {
		sweep(15092000+k);
		for (t=0; t<N; t++) {
			c_tilde[t][(k-k%D_BIN)/D_BIN]+=c(t)/(double)D_BIN; /*metodo intelligentissimo*/
		}
	}
	
	return c_tilde;

}

double* Cbar(double** c_tilde) {

	int t,i;
	double* c_bar=malloc(N*sizeof(double));
	
	

	for (t=0; t<N; t++) c_bar[t]=0;
	for (t=0; t<N; t++) {
		for (i=0; i<N_BIN; i++) {
			c_bar[t]+=c_tilde[t][i];
		}
		c_bar[t]/=(double)N_BIN;
	}

	return c_bar;
}




double ** Cluster (double ** c_tilde, double * c_bar) { /*Cluster(t,t_BIN)*/

	int t,i;
	double ** cluster=(double **) malloc(N*sizeof(double*));
	double *p=(double*)malloc(N_BIN*N*sizeof(double));
	for (t=0; t<N; t++) cluster[t]=p+t*N_BIN;
	/*fine dichiarazioni*/

	for (t=0; t<N; t++) {
		for (i=0; i<N_BIN; i++) {
			cluster[t][i]=c_bar[t]-(c_tilde[t][i]-c_bar[t])/(double)(N_BIN-1);
		}
	}

	return cluster;

}




double arcosh(double x) {
	return log(x+sqrt(x*x-1));
}

double* DeltaE(double * c_bar) { /*DeltaE(t)*/

	double *deltaE;
	int t;

	deltaE=malloc(N*sizeof(double));

	for (t=0; t<N; t++) 
		deltaE[t]=arcosh((c_bar[(t+1)%N]+c_bar[(t-1+N)%N])/(2*c_bar[t]));
	
	return deltaE;
}


double *Error_deltaE(double *deltaE, double **cluster) {

	int t,k;
	double * error_deltaE, **clusterE,*p;
	clusterE=(double **) malloc(N*sizeof(double*));
	p=(double*)malloc(N_BIN*N*sizeof(double));
	for (t=0; t<N; t++) clusterE[t]=p+t*N_BIN;
	error_deltaE=malloc(N*sizeof(double));
	/*fine dichiarazioni*/

	for (t=0; t<N; t++) {
		error_deltaE[t]=0;
		for (k=0; k<N_BIN; k++) {
			clusterE[t][k]=arcosh((cluster[(t+1)%N][k]+cluster[(t-1+N)%N][k])/(2*cluster[t][k]));
			error_deltaE[t]+=(double)(N_BIN-1)/(double)N_BIN*(clusterE[t][k]-deltaE[t])*(clusterE[t][k]-deltaE[t]);
		}
		error_deltaE[t]=sqrt(error_deltaE[t]);
	}

	return error_deltaE;


}





