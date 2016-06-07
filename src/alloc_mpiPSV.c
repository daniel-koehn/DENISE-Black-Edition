/*
 * Allocate memory for MPI variables (PSV problem) 
 *
 * Daniel Koehn
 * Kiel, 23/04/2016
 */

#include "fd.h"

void alloc_mpiPSV(struct mpiPSV *mpiPSV){

        /* global variables */
	extern int NX, NY, FDORDER;

	/* local variables */
	int nd, fdo3;

        nd = FDORDER/2 + 1;
        fdo3 = 2*nd;

	/* memory allocation for buffer arrays in which the wavefield 
   	information which is exchanged between neighbouring PEs is stored */
	(*mpiPSV).bufferlef_to_rig = matrix(1,NY,1,fdo3);
	(*mpiPSV).bufferrig_to_lef = matrix(1,NY,1,fdo3);
	(*mpiPSV).buffertop_to_bot = matrix(1,NX,1,fdo3);
	(*mpiPSV).bufferbot_to_top = matrix(1,NX,1,fdo3);

}



