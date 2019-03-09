/*------------------------------------------------------------------------
 * Output of RTM results for each shot (2D PSV problem)
 * 
 * Daniel Koehn
 * Kiel, 19/10/2016
 * ----------------------------------------------------------------------*/

#include "fd.h"

void RTM_PSV_out_shot(struct fwiPSV *fwiPSV, int ishot){

        /* global variables */
	extern int NX, NY, IDX, IDY;
	extern int POS[3], MYID_SHOT;
	extern char JACOBIAN[STRING_SIZE];
	extern MPI_Comm SHOT_COMM;
	
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

	MPI_Barrier(SHOT_COMM);

	/* merge P-image file */ 
	sprintf(jac1,"%s_P_image_shot_%i",JACOBIAN,ishot);
	if (MYID_SHOT==0) mergemod(jac1,3);

        /* clean up temporary files*/
        MPI_Barrier(SHOT_COMM);
        remove(jac);

	/* output of S-image */
        /* ----------------- */
	sprintf(jac,"%s_S_image_shot_%i.%i.%i",JACOBIAN,ishot,POS[1],POS[2]);
	FP=fopen(jac,"wb");

	for (i=1;i<=NX;i=i+IDX){
	   for (j=1;j<=NY;j=j+IDY){
                tmp = (*fwiPSV).waveconv_u_shot[j][i];
		fwrite(&tmp,sizeof(float),1,FP);
	   }
	}

	fclose(FP);

	MPI_Barrier(SHOT_COMM);

	/* merge S-image file */ 
	sprintf(jac1,"%s_S_image_shot_%i",JACOBIAN,ishot);
	if (MYID_SHOT==0) mergemod(jac1,3);

        /* clean up temporary files*/
        MPI_Barrier(SHOT_COMM);
        remove(jac);


}
