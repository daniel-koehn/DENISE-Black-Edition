/*
 * Allocate memory for material parameters (SH problem) 
 *
 * Daniel Koehn
 * Kiel, 03/12/2017
 */

#include "fd.h"

void alloc_matSH(struct matSH *matSH){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER;

	/* local variables */
	int nd;

        nd = FDORDER/2 + 1;	

	/* memory allocation for static (model) arrays */
	(*matSH).prho =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matSH).prhoi = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matSH).puip =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matSH).pujp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matSH).pu   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matSH).puipjp   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	/* memory allocation for static arrays for viscoelastic modeling */
	if (L>0){
		(*matSH).dip = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		(*matSH).d =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		(*matSH).e =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		(*matSH).ptaus =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
		(*matSH).ptausipjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
		(*matSH).fipjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
		(*matSH).f =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
		(*matSH).g =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
		(*matSH).peta =  vector(1,L);
		(*matSH).etaip =  vector(1,L);
		(*matSH).etajm =  vector(1,L);
		(*matSH).bip =  vector(1,L);
		(*matSH).bjm =  vector(1,L);
		(*matSH).cip =  vector(1,L);
		(*matSH).cjm =  vector(1,L);
	}
	
}



