/*
 *   Model homogeneous half space
 *   last update 11.04.02, T. Bohlen
 */

#include "fd.h"

void model_elastic(float  **  rho, float **  pi, float **  u){

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
	const float vp1=500.0, vs1=300.0, rho1=1800.0, h=15.0;
	
	/* parameters for layer 2 due to calculation of grad1, grad2 and grad3*/
	const float vp2=1200.0, vs2=700.0, rho2=2000.0;
	
	/*-----------------------------------------------------------------------*/

	y=h/DH;
	if(y==NYG) err(" \n y is equal NYG !! see src/model_grad.c  \n ");
	grad1=(vp2-vp1)/y;
	grad2=(vs2-vs1)/y;
	grad3=(rho2-rho1)/y;	
	
	
	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			
				if(j<=y){
				vp=vp1+(j*grad1);
				vs=vs1+(j*grad2);
				rhov=rho1+(j*grad3);
				}
				
				else{				
				vp=vp2;
				vs=vs2;
				rhov=rho2;
				}
				
				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;

					u[jj][ii]=vs;
					rho[jj][ii]=rhov;
					pi[jj][ii]=vp;
				}
			}
		}	

		
sprintf(modfile,"%s_rho_it_0.bin",INV_MODELFILE);
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
if (MYID==0) mergemod(modfile,3);
}



