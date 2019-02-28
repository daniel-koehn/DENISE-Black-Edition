/*
 * Allocate memory for FWI parameters (SH problem) 
 *
 * Daniel Koehn
 * Kiel, 13/12/2017
 */

#include "fd.h"

void alloc_fwiSH(struct fwiSH *fwiSH){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER, NTDTINV, NXNYI;

	/* local variables */
	int nd;

        nd = FDORDER/2 + 1;	

	/* elastic FWI parameters */
	(*fwiSH).prho_old =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).pu_old   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	
	(*fwiSH).Vs0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).Rho0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);	
	  
	(*fwiSH).gradg_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).gradp_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).waveconv_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).waveconv_rho_s = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).waveconv_rho_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	(*fwiSH).gradg_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).gradp_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).waveconv_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).waveconv_mu = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).waveconv_u_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	(*fwiSH).forward_prop_z =  vector(1,NXNYI*(NTDTINV));
	(*fwiSH).forward_prop_rho_z =  vector(1,NXNYI*(NTDTINV));
	(*fwiSH).forward_prop_sxz =  vector(1,NXNYI*(NTDTINV));
	(*fwiSH).forward_prop_syz =  vector(1,NXNYI*(NTDTINV));

	/* visco-elastic FWI parameters */
	/*if(L){

	    (*fwiPSV).Vp0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	    (*fwiPSV).ppi_old  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	    (*fwiPSV).waveconv = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	    (*fwiPSV).waveconv_lam = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	    (*fwiPSV).waveconv_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	    (*fwiPSV).gradg = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	    (*fwiPSV).gradp = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	}*/


}



