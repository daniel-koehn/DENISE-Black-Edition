/*
 * Rotate elastic tensor components (TTI PSV problem) 
 *
 * Daniel Koehn
 * Kiel, 06/02/2017
 */

#include "fd.h"

void rot_el_tensor_TTI(struct matTTI *matTTI){

        /* global variables */
	extern int NX, NY, L, FW, FDORDER;

	/* local variables */
	int i, j;
	float cos2, sin2, cos4, sin4, sincos, cos2sin2, cos2msin2;
	
	/* rotate elastic tensor components c11, c13, c33, c44 */
	for (j=1;j<=NY;j++){
	    for (i=1;i<=NX;i++){

		cos2 = cos((*matTTI).theta[j][i]) * cos((*matTTI).theta[j][i]);
		sin2 = sin((*matTTI).theta[j][i]) * sin((*matTTI).theta[j][i]);
		cos4 = cos2 * cos2;
		sin4 = sin2 * sin2;
		sincos = sin((*matTTI).theta[j][i]) * cos((*matTTI).theta[j][i]);
		cos2sin2 = cos2 * sin2;
		cos2msin2 = cos2 - sin2;

	        (*matTTI).d11[j][i] = cos2 * ((*matTTI).c11[j][i] * cos2 + (*matTTI).c13[j][i] * sin2) + sin2 * ((*matTTI).c13[j][i] * cos2 + (*matTTI).c33[j][i] * sin2) + 4.0 * (*matTTI).c44[j][i] * cos2sin2;
	        (*matTTI).d13[j][i] = (*matTTI).c11[j][i] * cos2sin2 + (*matTTI).c13[j][i] * sin4 + (*matTTI).c13[j][i] * cos4 + (*matTTI).c33[j][i] * cos2sin2 - 4.0 * (*matTTI).c44[j][i] * cos2sin2;
 	        (*matTTI).d15[j][i] = ((*matTTI).c13[j][i] - (*matTTI).c11[j][i]) * cos2 * sincos + ((*matTTI).c33[j][i] - (*matTTI).c13[j][i]) * sin2 * sincos + 2.0 * sincos * (*matTTI).c44[j][i] * cos2msin2;
	        (*matTTI).d33[j][i] = sin2 * ((*matTTI).c11[j][i] * sin2 + (*matTTI).c13[j][i] * cos2) + cos2 * ((*matTTI).c33[j][i] * cos2 + (*matTTI).c13[j][i] * sin2) + 4.0 * (*matTTI).c44[j][i] * cos2sin2;
 	        (*matTTI).d35[j][i] = ((*matTTI).c13[j][i] - (*matTTI).c11[j][i]) * sin2 * sincos + ((*matTTI).c33[j][i] - (*matTTI).c13[j][i]) * cos2 * sincos - 2.0 * sincos * (*matTTI).c44[j][i] * cos2msin2;
 	        (*matTTI).d55[j][i] = (*matTTI).c44[j][i] * (1.0 - 2.0 * sin2) * (1.0 - 2.0 * sin2) + ((*matTTI).c33[j][i] - (*matTTI).c13[j][i]) * cos2sin2 - ((*matTTI).c13[j][i] - (*matTTI).c11[j][i]) * cos2sin2;


	    }

	}	
	
}
