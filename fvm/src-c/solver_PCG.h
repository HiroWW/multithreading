#ifndef __H_SOLVER_PCG
#define __H_SOLVER_PCG
/*
extern int
solve_PCG(int N, int NL, int NU, int *INL, int **IAL, int *INU, int **IAU,
		double *D, double *B, double *X, double **AL, double **AU,
		double EPS, int *ITR, int *IER);
*/

extern int
solve_PCG(int N, int *indexLU, int *itemLU,
		double *D, double *B, double *X, double *AMAT, 
        	  double EPS, int *ITR, int *IER);

#endif /* __H_SOLVER_PCG */
