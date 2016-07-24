/*------------------------------------------------------------------------
 * Module for Limited Memory - Broyden-Fletcher-Goldfarb-Shanno (l-BFGS)
 * optimization (Nocedal and Wright, 2006)
 * 
 * Daniel Koehn
 * Kiel, 25.07.2016
 *
 * ----------------------------------------------------------------------*/

#include "fd.h"

void LBFGS(int iter, float * y_LBFGS, float * s_LBFGS, float * rho_LBFGS, float * alpha_LBFGS, float * q_LBFGS, float * r_LBFGS, float * beta_LBFGS, int LBFGS_pointer, int NLBFGS, int NLBFGS_vec){

	extern int MYID;
	
	int i, j, k, h, h1, h2;
	float gamma_LBFGS, sum_nom, sum_denom;
        float beta_LBFGS_1;
        int itershift;
	
        itershift = 1;

	/* update H^-1 * waveconv, using the l-BFGS method, if iter > 1 */
	/* ------------------------------------------------------------ */

	if(iter>1){
	     
	     /* calculate improved first guess Hessian gamma_LBFGS */
	     h1 = NLBFGS_vec*(LBFGS_pointer-1) + 1;
	     h2 = NLBFGS_vec*LBFGS_pointer; 
	     
	     sum_nom = dotp(y_LBFGS,s_LBFGS,h1,h2,0);
	     sum_denom = dotp(y_LBFGS,y_LBFGS,h1,h2,0);
	     gamma_LBFGS = sum_nom/sum_denom;
	     
	     /*printf("gamma_LBFGS = %e \n",gamma_LBFGS);*/
		 
	     /* update variable rho for all LBFGS-iterations and all parameter classes*/
	     for(k=1;k<=NLBFGS;k++){
		  
		h1 = NLBFGS_vec*(k-1) + 1;
		h2 = NLBFGS_vec*k;
		sum_nom = dotp(y_LBFGS,s_LBFGS,h1,h2,0); 
	
		if(fabs(sum_nom)>0.0){
		  rho_LBFGS[k] = 1.0/sum_nom;
		}
		else{
		  rho_LBFGS[k] = 0.0;
		} 
		  
		if(MYID==0){                                                
		printf("rho_LBFGS = %e of k = %d \n",rho_LBFGS[k],k);}
			                                               
	     }

	     /* update alpha_LBFGS and q_LBFGS */
	     for(k=NLBFGS;k>=1;k--){
		
	       h1 = NLBFGS_vec*(k-1) + 1;
	       h2 = NLBFGS_vec*k;
	       sum_nom = dotp(s_LBFGS,q_LBFGS,h1,h2,1);
	       alpha_LBFGS[k] = rho_LBFGS[k] * sum_nom;
	       
	       /* update q for all material parameters */
	       h = NLBFGS_vec*(k-1) + 1;
	       for (i=1;i<=NLBFGS_vec;i++){
		   q_LBFGS[i] = q_LBFGS[i] - alpha_LBFGS[k] * y_LBFGS[h];
		   h++;
	       }
	     }
		 
	       /* Multiply gradient with approximated Hessian */
	       for (i=1;i<=NLBFGS_vec;i++){
		   r_LBFGS[i] = gamma_LBFGS * q_LBFGS[i];
	       }

	     /* calculate H^-1 * waveconv[j][i] */
	     for(k=1;k<=NLBFGS;k++){
		
		h1 = NLBFGS_vec*(k-1) + 1;
		h2 = NLBFGS_vec*k;
		/* calculate beta_LBFGS*/   
		sum_nom = dotp(y_LBFGS,r_LBFGS,h1,h2,1);
		beta_LBFGS_1 = rho_LBFGS[k] * sum_nom;

		h = NLBFGS_vec*(k-1) + 1;
		for (i=1;i<=NLBFGS_vec;i++){
		   r_LBFGS[i] = r_LBFGS[i] + s_LBFGS[h]*(alpha_LBFGS[k]-beta_LBFGS_1);
		   h++;
		}
		 
	     }

	}
	
}
