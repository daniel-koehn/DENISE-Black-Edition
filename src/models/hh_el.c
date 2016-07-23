/*
 *   Model layer over half space (homogeneous)
 *   update 11.04.02, T. Bohlen
 *   last update 11. Okt. 2011 M. Schaefer
 */

#include "fd.h"

void model_elastic(float  **  rho, float **  pi, float **  u){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	extern char  MFILE[STRING_SIZE], INV_MODELFILE[STRING_SIZE];	
	extern float DH;
		/* local variables */
	float muv, piv, y, vp, vs, rhov;
	int i, j, ii, jj;
	 char modfile[STRING_SIZE];
	
	/* parameters for layer 1 */
	const float vp1=500.0, vs1=300.0, rho1=1800.0, h=5.0;
	
	/* parameters for layer 2 */
	const float vp2=500.0, vs2=300.0, rho2=1800.0;
	
	/*-----------------------------------------------------------------------*/



		

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			
				y=(float)j*DH;
			
				if (y<=h){
				 vp=vp1; vs=vs1; rhov=rho1; }

				
				else{
 				 vp=vp2; vs=vs2; rhov=rho2; }
				
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



