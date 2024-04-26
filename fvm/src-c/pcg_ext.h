#ifndef __H_PCG_EXT
#define __H_PCG_EXT

	extern int N2;
	extern int NLUmax, NLU;
	extern int METHOD, ORDER_METHOD;

	extern double EPSICCG;

        extern double *D, *PHI, *BFORCE;
        //extern double **AL, **AU;
        extern double *AMAT;

        extern int *INLU, *indexLU, *itemLU;
        extern int NPLU;

#endif /* __H_PCG_EXT */
