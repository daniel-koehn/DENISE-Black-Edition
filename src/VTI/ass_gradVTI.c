/*
 * Assemble RTM image for each shot (PSV VTI problem) 
 *
 * Daniel Koehn
 * Kiel, 15/02/2017
 */

#include "fd.h"

void ass_gradVTI(struct fwiPSV *fwiPSV, struct matVTI *matVTI, int iter){

        /* global variables */
	extern int NX, NY, IDX, IDY, INVMAT1;
        extern int GRAD_FORM;
        extern int INV_VP_ITER, INV_VS_ITER, INV_RHO_ITER;
	extern float DT;

	/* local variables */
	int i, j;

	/* FWI gradient calculation for VTI media not implemented yet */
	/* ---------------------------------------------------------- */	

  	/*for (i=1;i<=NX;i=i+IDX){
     	    for (j=1;j<=NY;j=j+IDY){
			 	
                if(iter<INV_VP_ITER){
           	     (*fwiPSV).waveconv_shot[j][i] = 0.0;
        	}
	                                                                       
      	    }
   	}*/
	
}



