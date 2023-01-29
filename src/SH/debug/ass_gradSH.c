/*
 * Assemble gradients for each shot (SH problem) 
 *
 * Daniel Koehn
 * Kiel, 13/12/2017
 */

#include "fd.h"

void ass_gradSH(struct fwiSH *fwiSH, struct matSH *matSH, int iter){

        /* global variables */
	extern int NX, NY, IDX, IDY, INVMAT1;
        extern int GRAD_FORM;
        extern int INV_VP_ITER, INV_VS_ITER, INV_RHO_ITER;
	extern float DT, C_vs, C_rho;

	/* local variables */
	int i, j;
        float muss, lamss;	

	/* calculate gradient for mu, Vs or Zs */
	/* ----------------------------------- */

	for (i=1;i<=NX;i=i+IDX){
   		for (j=1;j<=NY;j=j+IDY){
		 
      			/* calculate mu gradient */ 
      			(*fwiSH).waveconv_mu[j][i] = DT * (*fwiSH).waveconv_u_shot[j][i];
		 
      			if(INVMAT1==1){		
         			/* calculate Vs gradient */
				(*fwiSH).waveconv_u_shot[j][i] = 2.0 * (*matSH).prho[j][i] * (*matSH).pu[j][i] * (*fwiSH).waveconv_mu[j][i];
      			}
		 
      			if(INVMAT1==2){
        			/* calculate Zs gradient */
        			(*fwiSH).waveconv_u_shot[j][i] = 2.0 * (*matSH).pu[j][i] * (*fwiSH).waveconv_mu[j][i];
      			}
		 
      			if(INVMAT1==3){
        			/* calculate u gradient */
        			(*fwiSH).waveconv_u_shot[j][i] = (*fwiSH).waveconv_mu[j][i];
      			}

      			if(iter<INV_VS_ITER){
         			(*fwiSH).waveconv_u_shot[j][i] = 0.0;
      			}
		                                                                       
   		}
	}

	/* calculate gradient for density */
	/* ------------------------------ */
	for (i=1;i<=NX;i=i+IDX){
    		for (j=1;j<=NY;j=j+IDY){

       			/* calculate density gradient rho' */
       			(*fwiSH).waveconv_rho_s[j][i]= DT * (*fwiSH).waveconv_rho_shot[j][i];
				 
       			if(INVMAT1==1){
          			/* calculate density gradient */
          			(*fwiSH).waveconv_rho_shot[j][i] = ((*matSH).pu[j][i] * (*matSH).pu[j][i] * (*fwiSH).waveconv_mu[j][i]) + (*fwiSH).waveconv_rho_s[j][i];				 
       			}
		 
       			if(INVMAT1==3){
          			/* calculate density gradient */
          			(*fwiSH).waveconv_rho_shot[j][i] = (*fwiSH).waveconv_rho_s[j][i];
       			}

       			if(iter<INV_RHO_ITER){
          			(*fwiSH).waveconv_rho_shot[j][i] = 0.0;
       			}
	 
    		}
	}
	
}
