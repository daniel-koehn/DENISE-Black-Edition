/*
 * Allocate memory for FWI parameters (PSV problem) 
 *
 * Daniel Koehn
 * Kiel, 21/01/2016
 */

#include "fd.h"

void alloc_fwiPSV(struct fwiPSV *fwiPSV){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER, NTDTINV, NXNYI;

	/* local variables */
	int nd;

        nd = FDORDER/2 + 1;	

	(*fwiPSV).prho_old =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiPSV).ppi_old  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiPSV).pu_old   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	(*fwiPSV).Vp0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiPSV).Vs0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiPSV).Rho0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	(*fwiPSV).waveconv = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiPSV).waveconv_lam = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiPSV).waveconv_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	(*fwiPSV).gradg = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiPSV).gradp = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	  
	(*fwiPSV).gradg_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiPSV).gradp_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiPSV).waveconv_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiPSV).waveconv_rho_s = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiPSV).waveconv_rho_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	(*fwiPSV).gradg_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiPSV).gradp_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiPSV).waveconv_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiPSV).waveconv_mu = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiPSV).waveconv_u_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	(*fwiPSV).forward_prop_x =  vector(1,NXNYI*(NTDTINV));
	(*fwiPSV).forward_prop_y =  vector(1,NXNYI*(NTDTINV));
	(*fwiPSV).forward_prop_rho_x =  vector(1,NXNYI*(NTDTINV));
	(*fwiPSV).forward_prop_rho_y =  vector(1,NXNYI*(NTDTINV));
	(*fwiPSV).forward_prop_u =  vector(1,NXNYI*(NTDTINV));

}



