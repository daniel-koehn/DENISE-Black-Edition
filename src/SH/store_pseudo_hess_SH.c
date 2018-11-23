/*------------------------------------------------------------------------
 * Output of Pseudo-Hessian approximations
 * 
 * Daniel Koehn
 * Kiel, 27.10.2018
 *
 * ----------------------------------------------------------------------*/

#include "fd.h"

void store_pseudo_hess_SH(struct fwiSH *fwiSH){

	extern int NX, NY, IDX, IDY;
	extern int POS[3], MYID;
	extern char JACOBIAN[STRING_SIZE];
	
	char jac[225];
	int i, j;
	float tmp;
	FILE *FP1;

	/* ---------------------------------- */
	/* Write Pseudo-Hessian to hard drive */
	/* ---------------------------------- */

	/* Save main diagonal Hessian rho-rho */
	/* ---------------------------------- */
        sprintf(jac,"%s_hess_rho-rho.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP1=fopen(jac,"wb");
	
	for (i=1;i<=NX;i=i+IDX){   
           for (j=1;j<=NY;j=j+IDY){
		 tmp = (*fwiSH).hess_rho2p[j][i];
                 fwrite(&tmp,sizeof(float),1,FP1);
	   }
        }
        
	fclose(FP1);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* merge gradient file */ 
	sprintf(jac,"%s_hess_rho-rho",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
	/* Save main diagonal Hessian vs-vs */
	/* ---------------------------------- */
        sprintf(jac,"%s_hess_vs-vs.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP1=fopen(jac,"wb");
	
	for (i=1;i<=NX;i=i+IDX){   
           for (j=1;j<=NY;j=j+IDY){
		 tmp = (*fwiSH).hess_vs2[j][i];
                 fwrite(&tmp,sizeof(float),1,FP1);
	   }
        }
        
	fclose(FP1);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* merge gradient file */ 
	sprintf(jac,"%s_hess_vs-vs",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
	/* Save main diagonal Hessian ts-ts   */
	/* -------------------------------- */
        sprintf(jac,"%s_hess_ts-ts.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP1=fopen(jac,"wb");
	
	for (i=1;i<=NX;i=i+IDX){   
           for (j=1;j<=NY;j=j+IDY){
		 tmp = (*fwiSH).hess_ts2[j][i];
                 fwrite(&tmp,sizeof(float),1,FP1);
	   }
        }
        
	fclose(FP1);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* merge gradient file */ 
	sprintf(jac,"%s_hess_ts-ts",JACOBIAN);
	if (MYID==0) mergemod(jac,3);	
	
}
