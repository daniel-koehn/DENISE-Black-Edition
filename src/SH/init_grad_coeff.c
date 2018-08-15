/*
 * Calculate FWI gradient weighting coefficients (SH problem) 
 *
 * Daniel Koehn
 * Kiel, 14/08/2018
 */

#include "fd.h"

void init_grad_coeff(struct fwiSH *fwiSH, struct matSH *matSH){

        /* global variables */
	extern int NX, NY, L, ALPHA_VISC;

	/* local variables */	
	int i, j;

	for (j=1;j<=NY;j++){
	    for (i=1;i<=NX;i++){

		(*fwiSH).c1mu[j][i] = (1.0/((*matSH).pu[j][i]*(*matSH).pu[j][i])) * (1 + ALPHA_VISC * (*matSH).ptaus[j][i]) / (1 + L * (*matSH).ptaus[j][i]);
		(*fwiSH).c4mu[j][i] = (1.0/((*matSH).pu[j][i]*(*matSH).pu[j][i])) * (1 + ALPHA_VISC * (*matSH).ptaus[j][i]) / (*matSH).ptaus[j][i];

		(*fwiSH).c1ts[j][i] = (1.0/((*matSH).pu[j][i])) * (L - ALPHA_VISC) / ((1 + L * (*matSH).ptaus[j][i]) * (1 + L * (*matSH).ptaus[j][i]));
		(*fwiSH).c4ts[j][i] =  1.0/((*matSH).pu[j][i] * (*matSH).ptaus[j][i] * (*matSH).ptaus[j][i]);

	    }	
	}



}



