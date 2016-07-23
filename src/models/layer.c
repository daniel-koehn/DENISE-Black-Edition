/*
 *   Layer over halfspace
 *   M.Schäfer May 2011
 */

#include "fd.h"

void model_elastic(float  **  rho, float **  pi, float **  u){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DH;
	extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	extern char  MFILE[STRING_SIZE];	
	extern char INV_MODELFILE[STRING_SIZE];
		/* local variables */
	float vp, vs, rhov, y;
	int i, j, ii, jj;
	 
	
	/* parameters for layer 1 */
	const float vp1=680.0, vs1=320.0, rho1=1700.0, h=3.0;
	
	/* parameters for layer 2 */
	const float vp2=1000.0, vs2=590.0, rho2=2000.0;
	
	
	char modfile[STRING_SIZE];
	
	
	/*-----------------------------------------------------------------------*/



		

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			
				
				/* calculate coordinate in m */
				y=(float)j*DH;
				
 				/* two layer case */
 	                        if (y<=h){
																	                                   vp=vp1; vs=vs1; rhov=rho1; }


	                        else{
	                              vp=vp2; vs=vs2; rhov=rho2; }
			
				/*muv=vs1*vs1*rho1;
				piv=vp1*vp1*rho1;
				*/
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



