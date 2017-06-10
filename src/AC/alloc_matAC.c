/*
 * Allocate memory for material parameters (AC problem) 
 *
 * Daniel Koehn
 * Kiel, 10/06/2017
 */

#include "fd.h"

void alloc_matAC(struct matAC *matAC){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER;

	/* local variables */
	int nd;

        nd = FDORDER/2 + 1;	

	/* memory allocation for static (model) arrays */
	(*matAC).prho =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matAC).prip =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matAC).prjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matAC).ppi  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	/* memory allocation for static arrays for viscoelastic modeling */
	if (L>0){
		(*matAC).dip = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		(*matAC).d =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		(*matAC).e =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
		(*matAC).ptaup =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
		(*matAC).fipjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
		(*matAC).f =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
		(*matAC).g =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
		(*matAC).peta =  vector(1,L);
		(*matAC).etaip =  vector(1,L);
		(*matAC).etajm =  vector(1,L);
		(*matAC).bip =  vector(1,L);
		(*matAC).bjm =  vector(1,L);
		(*matAC).cip =  vector(1,L);
		(*matAC).cjm =  vector(1,L);
	}
	
}



