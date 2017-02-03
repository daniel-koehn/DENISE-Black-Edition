/*
 * Allocate memory for material parameters (VTI problem) 
 *
 * Daniel Koehn
 * Kiel, 01/02/2017
 */

#include "fd.h"

void alloc_matVTI(struct matVTI *matVTI){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER;

	/* local variables */
	int nd;

        nd = FDORDER/2 + 1;	

	/* memory allocation for static (model) arrays */
	(*matVTI).prho =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matVTI).prip =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matVTI).prjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matVTI).c11  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matVTI).c13   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matVTI).c33   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matVTI).c44   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matVTI).c44h   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	
}



