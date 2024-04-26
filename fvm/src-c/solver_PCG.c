/*
 * solver_PCG
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "solver_PCG.h"

extern int
solve_PCG(int N, int *indexLU, int *itemLU,
		double *D, double *B, double *X, double *AMAT, 
		double EPS, int *ITR, int *IER)
{
	double **W;
	double VAL, BNRM2, WVAL, SW, RHO, BETA, RHO1, C1, DNRM2, ALPHA, ERR;
	double Stime, Etime;	
	int i, j, ic, ip, L, ip1;
	int R = 0;
	int Z = 1;
	int Q = 1;
	int P = 2;
	int DD = 3;

/*********
 * INIT. *
 *********/

	W = (double **)malloc(sizeof(double *)*4);
	if(W == NULL) {
		fprintf(stderr, "Error: %s\n", strerror(errno));
	return -1;
	}
	for(i=0; i<4; i++) {
		W[i] = (double *)malloc(sizeof(double)*N);
		if(W[i] == NULL) {
			fprintf(stderr, "Error: %s\n", strerror(errno));
			return -1;
		}
	}

	for(i=0; i<N; i++) {
		X[i] = 0.0;
		W[1][i] = 0.0;
		W[2][i] = 0.0;
		W[3][i] = 0.0;
	}

	  for(i=0; i<N; i++) {
	    W[DD][i] = 1.0 / D[i];
	  }

/**************************
 * {r0} = {b} - {A}{xini} *
 **************************/
	for(i=0; i<N; i++) {
		VAL = D[i] * X[i];
		for(j=indexLU[i]; j<indexLU[i+1]; j++) {
			VAL += AMAT[j] * X[itemLU[j]];
		}
		W[R][i] = B[i] - VAL;
	}

	BNRM2 = 0.0;
	for(i=0; i<N; i++) {
	  BNRM2 += B[i]*B[i];	  
	}

/************************************************************** ITERATION */
	*ITR = N;

	for(L=0; L<(*ITR); L++) {

/*******************
 * {z} = [Minv]{r} *
 *******************/
	  for(i=0; i<N; i++) {
	    W[Z][i] = W[R][i]*W[DD][i];
	  }

/****************
 * RHO = {r}{z} *
 ****************/
		RHO = 0.0;
		for(i=0; i<N; i++) {
		  RHO += W[R][i] * W[Z][i];
		}

/********************************
 * {p}  = {z} if      ITER=0    *
 * BETA = RHO / RHO1  otherwise *
 ********************************/
		if(L == 0) {
		  for(i=0; i<N; i++) {
			W[P][i] = W[Z][i];
		  }
		} else {
		  BETA = RHO / RHO1;
		  for(i=0; i<N; i++) {
			W[P][i] = W[Z][i] + BETA * W[P][i];
		}
		}

/****************
 * {q} = [A]{p} *
 ****************/
		for(i=0; i<N; i++) {
		  VAL = D[i] * W[P][i];
		  for(j=indexLU[i]; j<indexLU[i+1]; j++) {
			VAL += AMAT[j] * W[P][itemLU[j]];
		  }
		  W[Q][i] = VAL;
		}

/************************
 * ALPHA = RHO / {p}{q} *
 ************************/
		C1 = 0.0;
		for(i=0; i<N; i++) {
			C1 += W[P][i] * W[Q][i];
		}

		ALPHA = RHO / C1;

/***************************
 * {x} = {x} + ALPHA * {p} *
 * {r} = {r} - ALPHA * {q} *
 ***************************/
		for(i=0; i<N; i++) {
			X[i]    += ALPHA * W[P][i];
			W[R][i] -= ALPHA * W[Q][i];
		}

		DNRM2 = 0.0;
		for(i=0; i<N; i++) {
			DNRM2 += pow(W[R][i], 2.0);
		}

		ERR = sqrt(DNRM2/BNRM2);
                if( (L+1)%100 ==1) {
                        fprintf(stderr, "%5d%16.6e\n", L+1, ERR);
                }

		if(ERR < EPS) {
			*IER = 0;
			goto N900;
		} else {
			RHO1 = RHO;
		}
	}
	*IER = 1;

N900:
	fprintf(stdout, "%5d%16.6e\n", L+1, ERR);
	*ITR = L;

	free(W);

	return 0;
}