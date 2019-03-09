/*
 * Calculate FWI gradient weighting coefficients (SH problem) 
 *
 * Daniel Koehn
 * Kiel, 14/08/2018
 */

#include "fd.h"

void init_grad_coeff(struct fwiSH *fwiSH, struct matSH *matSH){

        /* global variables */
	extern int NX, NY, L, MYID;
	extern float ALPHA_VISC, C_vs, C_rho, C_taus;

	/* local variables */	
	int i, j;
	float C_mu, C_taus1;

	/* Calculate shear modulus */
	// C_mu = C_rho * C_vs * C_vs;

	for (j=1;j<=NY;j++){
	    for (i=1;i<=NX;i++){

		C_mu = (*matSH).pu[j][i] * (*matSH).pu[j][i] * (*matSH).prho[j][i];
		C_taus1 = (*matSH).ptaus[j][i];

                (*fwiSH).c1mu[j][i] = (1.0/(C_mu * C_mu)) * (1 + ALPHA_VISC * C_taus1) / (1 + L * C_taus1);
		(*fwiSH).c4mu[j][i] = (1.0/(C_mu * C_mu)) * (1 + ALPHA_VISC * C_taus1) / C_taus1;

		(*fwiSH).c1ts[j][i] = (1.0/C_mu) * (L - ALPHA_VISC) / ((1 + L * C_taus1) * (1 + L * C_taus1));
		(*fwiSH).c4ts[j][i] =  1.0/(C_mu * C_taus1 * C_taus1);

	    }	
	}
	
	if(MYID==0){
	    printf("\n");
	    printf("Compute gradient coefficients c1mu, c1ts, c4mu, c4ts\n");
	    printf("----------------------------------------------------\n");
	    printf("c1mu = %e \t c4mu = %e \n",(*fwiSH).c1mu[1][1],(*fwiSH).c4mu[1][1]);
	    printf("c1ts = %e \t c4ts = %e \n",(*fwiSH).c1ts[1][1],(*fwiSH).c4ts[1][1]);
	    printf("\n");
	}

}



