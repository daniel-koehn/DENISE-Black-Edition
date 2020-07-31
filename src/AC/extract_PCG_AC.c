/*------------------------------------------------------------------------
 * Module for extraction of Preconditioned Conjugate Gradient Method (PCG)
 * for the material parameters vp, rho and lambda, rho respectively
 * 
 * Daniel Koehn
 * Kiel, 12.06.2017
 * ----------------------------------------------------------------------*/

#include "fd.h"

void extract_PCG_AC(float * PCG_old, float ** waveconv, float ** waveconv_rho){

	extern int NX, NY, IDX, IDY, POS[3], MYID;
	extern char JACOBIAN[STRING_SIZE];
	
	int i, j, h;
	char jac[STRING_SIZE];
	FILE *FP3;
	
	
	/* ============================================================================================================================================================== */
	/* ===================================================== GRADIENT VP/ZP/lambda ================================================================================== */
	/* ============================================================================================================================================================== */

        h=1;
	/* store gradient */
	for (i=1;i<=NX;i=i+IDX){
	   for (j=1;j<=NY;j=j+IDY){

		 waveconv[j][i] = PCG_old[h];

                 h++;
	   }
	}
	
	/* save gradient */
	sprintf(jac,"%s_p.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                	fwrite(&waveconv[j][i],sizeof(float),1,FP3);
            }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
	/* clean up temporary files from each MPI process */
        MPI_Barrier(MPI_COMM_WORLD);
        sprintf(jac,"%s_p.%i.%i",JACOBIAN,POS[1],POS[2]);
        remove(jac);

	/* ============================================================================================================================================================== */
	/* ===================================================== GRADIENT rho =========================================================================================== */
	/* ============================================================================================================================================================== */

	/* store gradient */
	for (i=1;i<=NX;i=i+IDX){
	   for (j=1;j<=NY;j=j+IDY){

		 waveconv_rho[j][i] = PCG_old[h];

                 h++;
	   }
	}
	
	/* save old gradient */
	sprintf(jac,"%s_p_rho.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                	fwrite(&waveconv_rho[j][i],sizeof(float),1,FP3);
            }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_rho",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
	/* clean up temporary files from each MPI process */
        MPI_Barrier(MPI_COMM_WORLD);
        sprintf(jac,"%s_p_rho.%i.%i",JACOBIAN,POS[1],POS[2]);
        remove(jac);


}
