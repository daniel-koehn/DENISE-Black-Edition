/*------------------------------------------------------------------------
 * Output of RTM results for each shot (2D SH problem)
 * 
 * Daniel Koehn
 * Kiel, 13/12/2017
 * ----------------------------------------------------------------------*/

#include "fd.h"

void RTM_SH_out_shot(struct fwiSH *fwiSH, int ishot){

        /* global variables */
	extern int NX, NY, IDX, IDY;
	extern int POS[3], MYID;
	extern char JACOBIAN[STRING_SIZE];
	
        /* local variables */
	char jac[STRING_SIZE], jac1[STRING_SIZE];
	int i, j;
        float tmp;
	FILE *FP;	

	/* output of S-image */
        /* ----------------- */
	sprintf(jac,"%s_S_image_shot_%i.%i.%i",JACOBIAN,ishot,POS[1],POS[2]);
	FP=fopen(jac,"wb");

	for (i=1;i<=NX;i=i+IDX){
	   for (j=1;j<=NY;j=j+IDY){
                tmp = (*fwiSH).waveconv_u_shot[j][i];
		fwrite(&tmp,sizeof(float),1,FP);
	   }
	}

	fclose(FP);

	MPI_Barrier(MPI_COMM_WORLD);

	/* merge S-image file */ 
	sprintf(jac1,"%s_S_image_shot_%i",JACOBIAN,ishot);
	if (MYID==0) mergemod(jac1,3);

        /* clean up temporary files*/
        MPI_Barrier(MPI_COMM_WORLD);
        remove(jac);

}
