/*------------------------------------------------------------------------
 *   zero LBFGS-vectors
 *  
 *  
 *   last update 05/08/2013, D. Koehn
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void zero_LBFGS(int NLBFGS, int NLBFGS_vec, float * y_LBFGS, float * s_LBFGS, float * q_LBFGS, float * r_LBFGS, 
                 float * alpha_LBFGS, float * beta_LBFGS, float * rho_LBFGS){


	register int i;

	for (i=1;i<=(NLBFGS_vec*NLBFGS);i++){
            y_LBFGS[i]=0.0;
            s_LBFGS[i]=0.0;
	}

	for (i=1;i<=(NLBFGS_vec);i++){
            q_LBFGS[i]=0.0;
            r_LBFGS[i]=0.0;
	}
	
        for (i=1;i<=(NLBFGS);i++){
            alpha_LBFGS[i]=0.0;
            beta_LBFGS[i]=0.0;
            rho_LBFGS[i]=0.0;
	}	
	
}
