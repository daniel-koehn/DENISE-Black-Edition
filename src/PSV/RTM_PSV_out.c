/*------------------------------------------------------------------------
 * Output of RTM results (2D PSV problem)
 * 
 * Daniel Koehn
 * Kiel, 06/06/2016
 * ----------------------------------------------------------------------*/

#include "fd.h"

void RTM_PSV_out(struct fwiPSV *fwiPSV){

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
	sprintf(jac,"%s_P_image.%i%i",JACOBIAN,POS[1],POS[2]);
	FP=fopen(jac,"wb");

	for (i=1;i<=NX;i=i+IDX){
	   for (j=1;j<=NY;j=j+IDY){
                tmp = (*fwiPSV).waveconv[j][i];
		fwrite(&tmp,sizeof(float),1,FP);
	   }
	}

	fclose(FP);

	MPI_Barrier(MPI_COMM_WORLD);

	/* merge P-image file */ 
	sprintf(jac1,"%s_P_image",JACOBIAN);
	if (MYID==0) mergemod(jac1,3);

        /* clean up temporary files*/
        MPI_Barrier(MPI_COMM_WORLD);
        remove(jac);

	/* output of S-image */
        /* ----------------- */
	sprintf(jac,"%s_S_image.%i%i",JACOBIAN,POS[1],POS[2]);
	FP=fopen(jac,"wb");

	for (i=1;i<=NX;i=i+IDX){
	   for (j=1;j<=NY;j=j+IDY){
                tmp = (*fwiPSV).waveconv_u[j][i];
		fwrite(&tmp,sizeof(float),1,FP);
	   }
	}

	fclose(FP);

	MPI_Barrier(MPI_COMM_WORLD);

	/* merge S-image file */ 
	sprintf(jac1,"%s_S_image",JACOBIAN);
	if (MYID==0) mergemod(jac1,3);

        /* clean up temporary files*/
        MPI_Barrier(MPI_COMM_WORLD);
        remove(jac);


}
