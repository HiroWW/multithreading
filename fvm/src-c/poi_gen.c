/*
 * POI_GEN
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "struct_ext.h"
#include "pcg_ext.h"
#include "poi_gen.h"
#include "allocate.h"

extern int
POI_GEN(void)
{
	int nn;
	int ic0, icN1, icN2, icN3, icN4, icN5, icN6;
	int i, j, k, ib, ic, ip, icel, icou, icol, icouG;
	int ii, jj, kk, nn1, num, nr, j0, j1;
	double coef, VOL0, S1t, E1t;
        int isLU;


	NLU= 6;

	BFORCE = (double *)allocate_vector(sizeof(double),ICELTOT);
	D      = (double *)allocate_vector(sizeof(double),ICELTOT);
	PHI    = (double *)allocate_vector(sizeof(double),ICELTOT);
        INLU = (int *)allocate_vector(sizeof(int),ICELTOT);

        indexLU= (int *)allocate_vector(sizeof(int),ICELTOT+1);

        for (i = 0; i <ICELTOT ; i++) {
		BFORCE[i]=0.0;
		D[i]  =0.0;
		PHI[i]=0.0;
        	INLU[i] = 0;
        }

        for (i = 0; i <=ICELTOT ; i++) {
		indexLU[i] = 0;
	}


/*********************************
 * INTERIOR & NEUMANN boundary's *
 *********************************/

	for(icel=0; icel<ICELTOT; icel++) {
		icN1 = NEIBcell[icel][0];
		icN2 = NEIBcell[icel][1];
		icN3 = NEIBcell[icel][2];
		icN4 = NEIBcell[icel][3];
		icN5 = NEIBcell[icel][4];
		icN6 = NEIBcell[icel][5];

		if(icN5 != 0) {
 		   INLU[icel]= INLU[icel] + 1;		  
		}

		if(icN3 != 0) {
 		   INLU[icel]= INLU[icel] + 1;
		}

		if(icN1 != 0) {
 		   INLU[icel]= INLU[icel] + 1;		  
		}

		if(icN2 != 0) {
 		   INLU[icel]= INLU[icel] + 1;		  
		}

		if(icN4 != 0) {
 		   INLU[icel]= INLU[icel] + 1;		  
		}

		if(icN6 != 0) {
 		   INLU[icel]= INLU[icel] + 1;		  
		}
	}


/********************************************
* 1D ordering: indexL, indexU, itemL, itemU *
*********************************************/

        for(i=0; i<ICELTOT; i++){
                indexLU[i+1]=indexLU[i]+INLU[i];
        }
        NPLU= indexLU[ICELTOT];

        itemLU= (int *)allocate_vector(sizeof(int),NPLU);
        AMAT  = (double *)allocate_vector(sizeof(double),NPLU);

	for(i=0; i<ICELTOT; i++) {
	  for(j=indexLU[i]; j<indexLU[i+1]; j++) {
	    itemLU[j]=0;
	      AMAT[j]=0.0;
	  }
	}

        free(INLU);

/*************************************
 * INTERIOR & NEUMANN BOUNDARY CELLs *
 *************************************/
	for(icel=0; icel<ICELTOT; icel++) {
			icN1 = NEIBcell[icel][0];
			icN2 = NEIBcell[icel][1];
			icN3 = NEIBcell[icel][2];
			icN4 = NEIBcell[icel][3];
			icN5 = NEIBcell[icel][4];
			icN6 = NEIBcell[icel][5];

	                VOL0 = VOLCEL[icel];

       		        isLU=indexLU[icel];

                        icou= 0;
			if(icN5 != 0) {
				coef = RDZ * ZAREA;
				D[icel] -= coef;
				itemLU[icou+isLU]= icN5-1;
				AMAT  [icou+isLU]= coef;
				icou= icou + 1;
			}

			if(icN3 != 0) {
				coef = RDZ * YAREA;
				D[icel] -= coef;
				itemLU[icou+isLU]= icN3-1;
				AMAT  [icou+isLU]= coef;
				icou= icou + 1;
			}

			if(icN1 != 0) {
				coef = RDZ * XAREA;
				D[icel] -= coef;
				itemLU[icou+isLU]= icN1-1;
				AMAT  [icou+isLU]= coef;
				icou= icou + 1;
			}

			if(icN2 != 0) {
				coef = RDZ * XAREA;
				D[icel] -= coef;
				itemLU[icou+isLU]= icN2-1;
				AMAT  [icou+isLU]= coef;
				icou= icou + 1;
			}

			if(icN4 != 0) {
				coef = RDZ * YAREA;
				D[icel] -= coef;
				itemLU[icou+isLU]= icN4-1;
				AMAT  [icou+isLU]= coef;
				icou= icou + 1;
			}

			if(icN6 != 0) {
			        coef = RDZ * ZAREA;
				D[icel] -= coef;
				itemLU[icou+isLU]= icN6-1;
				AMAT  [icou+isLU]= coef;
				icou= icou + 1;
			}

			ii = XYZ[icel][0];
			jj = XYZ[icel][1];
			kk = XYZ[icel][2];

			BFORCE[icel] = - (double)(ii + jj + kk) * VOLCEL[icel];
		}

/****************************
 * DIRICHLET BOUNDARY CELLs *
 ****************************/
/* TOP SURFACE */
	for(ib=0; ib<ZmaxCELtot; ib++) {
		icel  = ZmaxCEL[ib];
		coef = 2.0 * RDZ * ZAREA;
		D[icel-1] -= coef;
	}


	return 0;
}
