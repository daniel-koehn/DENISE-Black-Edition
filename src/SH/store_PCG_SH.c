/*------------------------------------------------------------------------
 * Module for storage of Preconditioned Conjugate Gradient Method (PCG)
 * for the material parameters vs, rho and mu, rho respectively
 * 
 * Daniel Koehn
 * Kiel, 22.12.2017
 * ----------------------------------------------------------------------*/

#include "fd.h"

void store_PCG_SH(float * PCG_old, float ** waveconv_u, float ** waveconv_rho, float ** waveconv_ts){

	extern int NX, NY, IDX, IDY;
	extern int POS[3], MYID;
	extern char JACOBIAN[STRING_SIZE];
	
	char jac[225], jac1[225];
	int i, j, h;
	FILE *FP3;
	
	/* ============================================================================================================================================================== */
	/* ===================================================== GRADIENT VS/ZS/mu ====================================================================================== */
	/* ============================================================================================================================================================== */
	
	h=1;
	/* store gradient */
	for (i=1;i<=NX;i=i+IDX){
	   for (j=1;j<=NY;j=j+IDY){

		 PCG_old[h] = waveconv_u[j][i];

                 h++;
	   }
	}

	/* ============================================================================================================================================================== */
	/* ===================================================== GRADIENT rho =========================================================================================== */
	/* ============================================================================================================================================================== */

	/* store gradient */
	for (i=1;i<=NX;i=i+IDX){
	   for (j=1;j<=NY;j=j+IDY){

		 PCG_old[h] = waveconv_rho[j][i];

                 h++;
	   }
	}

	/* ============================================================================================================================================================== */
	/* ===================================================== GRADIENT rho =========================================================================================== */
	/* ============================================================================================================================================================== */

	/* store gradient */
	for (i=1;i<=NX;i=i+IDX){
	   for (j=1;j<=NY;j=j+IDY){

		 PCG_old[h] = waveconv_ts[j][i];

                 h++;
	   }
	}	

	/* save Vs gradient */
	/* ---------------- */
        sprintf(jac1,"%s_p_u.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac1,"wb");
	
	for (i=1;i<=NX;i=i+IDX){   
           for (j=1;j<=NY;j=j+IDY){
                 fwrite(&waveconv_u[j][i],sizeof(float),1,FP3);
	   }
        }
        
	fclose(FP3);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* merge gradient file */ 
	sprintf(jac,"%s_p_u.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	//remove(jac1);

	/* save density gradient */
	/* --------------------- */
        sprintf(jac1,"%s_p_rho.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac1,"wb");
	
	for (i=1;i<=NX;i=i+IDX){   
           for (j=1;j<=NY;j=j+IDY){
                 fwrite(&waveconv_rho[j][i],sizeof(float),1,FP3);
	   }
        }
        
	fclose(FP3);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* merge gradient file */ 
	sprintf(jac,"%s_p_rho.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	//remove(jac1);

	/* save Taus gradient */
	/* ------------------ */
        sprintf(jac1,"%s_p_taus.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac1,"wb");
	
	for (i=1;i<=NX;i=i+IDX){   
           for (j=1;j<=NY;j=j+IDY){
                 fwrite(&waveconv_ts[j][i],sizeof(float),1,FP3);
	   }
        }
        
	fclose(FP3);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* merge gradient file */ 
	sprintf(jac,"%s_p_taus.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	//remove(jac1);
	
}
