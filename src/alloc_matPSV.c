/*
 * Allocate memory for material parameters (PSV problem) 
 *
 * Daniel Koehn
 * Kiel, 21/01/2016
 */

#include "fd.h"

void alloc_matPSV(struct matPSV *matPSV){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER;

	/* local variables */
	int nd;

        nd = FDORDER/2 + 1;	

	/* memory allocation for static (model) arrays */
	(*matPSV).prho =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matPSV).prip =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matPSV).prjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matPSV).ppi  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matPSV).pu   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matPSV).puipjp   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	/* memory allocation for static arrays for viscoelastic modeling */
	if (L>0){
		(*matPSV).dip = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		(*matPSV).d =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		(*matPSV).e =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		(*matPSV).ptaus =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
		(*matPSV).ptausipjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
		(*matPSV).ptaup =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
		(*matPSV).fipjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
		(*matPSV).f =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
		(*matPSV).g =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
		(*matPSV).peta =  vector(1,L);
		(*matPSV).etaip =  vector(1,L);
		(*matPSV).etajm =  vector(1,L);
		(*matPSV).bip =  vector(1,L);
		(*matPSV).bjm =  vector(1,L);
		(*matPSV).cip =  vector(1,L);
		(*matPSV).cjm =  vector(1,L);
	}
	
}



