/*
 * Allocate memory for material parameters (TTI problem) 
 *
 * Daniel Koehn
 * Kiel, 04/02/2017
 */

#include "fd.h"

void alloc_matTTI(struct matTTI *matTTI){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER;

	/* local variables */
	int nd;

        nd = FDORDER/2 + 1;	

	/* memory allocation for static (model) arrays */
	(*matTTI).prho =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matTTI).prip =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matTTI).prjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matTTI).c11  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matTTI).c13   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matTTI).c33   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matTTI).c44   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	(*matTTI).d11  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matTTI).d13  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matTTI).d15  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matTTI).d33  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matTTI).d35  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matTTI).d55  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	(*matTTI).d15h  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matTTI).d35h  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	(*matTTI).d55h  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);

	(*matTTI).theta  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	
}
