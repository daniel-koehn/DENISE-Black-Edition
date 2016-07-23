/*
 *   Model homogeneous half space
 *   last update 23.04.03, T. Bohlen
 */

#include "fd.h"

void model_elastic(float  **  rho, float **  pi, float **  u){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DH, REC_ARRAY_DEPTH;
	extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	extern char  MFILE[STRING_SIZE];	

		/* local variables */
	float rhov, muv, piv, vp, vs;
	float y, x, y0, ys;
	int i, j, ii, jj;
	 
	
	/* parameters for layer 1 */
	/*const float vp1=0.0, vs1=0.0001, rho1=1.25, h=10.0, xs=30.0, phi=PI*REC_ARRAY_DEPTH/180.0;*/

	const float vp1=1700.0, vs1=170.0, rho1=1700.0, h=10.0, xs=30.0, phi=PI*REC_ARRAY_DEPTH/180.0;

	/* parameters for layer 2 */
	/*const float vp2=5700.0, vs2=3400.0, rho2=2200.0;*/
	const float vp2=1700.0, vs2=170.0, rho2=1700.0;
	
	
	/*-----------------------------------------------------------------------*/



		

	y0=h-xs*tan(phi);
	
	/* loop over global grid */
		for (i=-1;i<=NXG+2;i++){
			for (j=-1;j<=NYG+2;j++){
			
				x=(float)i*DH;
				y=(float)j*DH;
				
				
				ys=y0+tan(phi)*x;
			
				if (y<ys){
				 vp=vp1; vs=vs1; rhov=rho1; }

				
				 else{
 				 vp=vp2; vs=vs2; rhov=rho2; }
                    
				
				muv=vs*vs*rhov;
				piv=vp*vp*rhov;

				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;

					
					u[jj][ii]=muv;
					rho[jj][ii]=rhov;
					pi[jj][ii]=piv;
				}
			}
		}	

		

	
	/* each PE writes his model to disk */

	writemod(MFILE,pi,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(MFILE,3);

}



