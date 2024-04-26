#ifndef __H_PCG
#define __H_PCG

        int N2;
        int NLUmax, NLU;
	int METHOD, ORDER_METHOD;

	double EPSICCG;

	double *D, *PHI, *BFORCE;
	double *AMAT;

        int *INLU, *indexLU, *itemLU;
	int NPLU;

#endif /* __H_PCG */
