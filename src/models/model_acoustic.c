/*
 *   Model homogeneous half space
 *   last update 11.06.2017, T. Bohlen
 */

#include "fd.h"

void model_acoustic(float  **  rho, float **  pi){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	extern char  MFILE[STRING_SIZE];	
	extern char INV_MODELFILE[STRING_SIZE];
	extern float DH;
		/* local variables */
	float vp, vs, rhov, grad1, grad2, grad3, y;
	int i, j, ii, jj;
	char modfile[STRING_SIZE]; 
	
	/* parameters for layer 1 */
	const float vp1=3000.0, rho1=1000.0;
	
	/*-----------------------------------------------------------------------*/
	
	
	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
							
				vp=vp1;				
				rhov=rho1;
					
				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
					
					rho[jj][ii]=rhov;
					pi[jj][ii]=vp;
				}
			}
		}	

		
/*sprintf(modfile,"%s_rho_it_0.bin",INV_MODELFILE);
writemod(modfile,rho,3);
MPI_Barrier(MPI_COMM_WORLD);
if (MYID==0) mergemod(modfile,3);

sprintf(modfile,"%s_vs_it_0.bin",INV_MODELFILE);
writemod(modfile,u,3);
MPI_Barrier(MPI_COMM_WORLD);
if (MYID==0) mergemod(modfile,3);

sprintf(modfile,"%s_vp_it_0.bin",INV_MODELFILE);
writemod(modfile,pi,3);
MPI_Barrier(MPI_COMM_WORLD);
if (MYID==0) mergemod(modfile,3);*/

}



