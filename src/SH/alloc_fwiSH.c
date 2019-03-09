/*
 * Allocate memory for FWI parameters (SH problem) 
 *
 * Daniel Koehn
 * Kiel, 13/12/2017
 */

#include "fd.h"

void alloc_fwiSH(struct fwiSH *fwiSH){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER, NTDTINV, NXNYI, EPRECOND;

	/* local variables */
	int nd;

        nd = FDORDER/2 + 1;	

	/* (visco)elastic FWI parameters */

	/* temporary arrays to store material parameters */
	(*fwiSH).prho_old =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).pu_old   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).ptaus_old   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	
	/* material parameters for initial model */
	(*fwiSH).Vs0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).Rho0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);	
	(*fwiSH).Taus0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);	

	/* (preconditioned) gradients / Hessian-gradient products */
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

	(*fwiSH).gradg_ts = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).gradp_ts = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).waveconv_ts = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).waveconv_ts_s = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).waveconv_ts_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	/* arrays to store forward wavefields */
	(*fwiSH).forward_prop_z =  vector(1,NXNYI*(NTDTINV));
	(*fwiSH).forward_prop_rho_z =  vector(1,NXNYI*(NTDTINV));
	(*fwiSH).forward_prop_sxz =  vector(1,NXNYI*(NTDTINV));
	(*fwiSH).forward_prop_syz =  vector(1,NXNYI*(NTDTINV));

	(*fwiSH).forward_prop_rxz =  matrix(1,NXNYI*(NTDTINV),1,L);
	(*fwiSH).forward_prop_ryz =  matrix(1,NXNYI*(NTDTINV),1,L);

	/* weighting constants for gradient computation */
	(*fwiSH).c1mu = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).c4mu = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).c1ts = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).c4ts = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	/* time integrated memory variables */
	(*fwiSH).Rxz = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	(*fwiSH).Ryz = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);

	/* time derivative of memory variables */
	(*fwiSH).rxzt = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	(*fwiSH).ryzt = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);

	/* tausl */
	(*fwiSH).tausl =  vector(1,L);

	/* Pseudo Hessian main- and non-diagonal approximations */
	(*fwiSH).hess_mu2 = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).hess_rho2 = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).hess_ts2 = matrix(-nd+1,NY+nd,-nd+1,NX+nd);	
	(*fwiSH).hess_vs2 = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).hess_rho2p = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	(*fwiSH).hess_muts = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).hess_murho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*fwiSH).hess_tsrho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);	

}



