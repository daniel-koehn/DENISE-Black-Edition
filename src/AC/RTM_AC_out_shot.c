/*------------------------------------------------------------------------
 * Output of RTM results for each shot (2D AC problem)
 * 
 * Daniel Koehn
 * Kiel, 14/07/2017
 * ----------------------------------------------------------------------*/

#include "fd.h"

void RTM_AC_out_shot(struct fwiPSV *fwiPSV, int ishot){

        /* global variables */
	extern int NX, NY, IDX, IDY;
	extern int POS[3], MYID;
	extern char JACOBIAN[STRING_SIZE];
	
        /* local variables */
	char jac[STRING_SIZE], jac1[STRING_SIZE];
	int i, j;
        float tmp;
	FILE *FP;
	
	/* output of P-image */
        /* ----------------- */
	sprintf(jac,"%s_P_image_shot_%i.%i.%i",JACOBIAN,ishot,POS[1],POS[2]);
	FP=fopen(jac,"wb");

	for (i=1;i<=NX;i=i+IDX){
	   for (j=1;j<=NY;j=j+IDY){
                tmp = (*fwiPSV).waveconv_shot[j][i];
		fwrite(&tmp,sizeof(float),1,FP);
	   }
	}

	fclose(FP);

	MPI_Barrier(MPI_COMM_WORLD);

	/* merge P-image file */ 
	sprintf(jac1,"%s_P_image_shot_%i",JACOBIAN,ishot);
	if (MYID==0) mergemod(jac1,3);

        /* clean up temporary files*/
        MPI_Barrier(MPI_COMM_WORLD);
        remove(jac);

}
