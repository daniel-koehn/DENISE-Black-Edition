/*
 *   Model with slow formation and fluid-filled borehole
 *   last update 11.12.07, O. Hellwig
 */

#include "fd.h"

void model_elastic(float  **  rho, float **  pi, float **  u){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, DH;
	extern int   NX, NY, NXG, NYG,  POS[3], MYID;
	extern char  MFILE[STRING_SIZE];	

	/* local variables */
	float rhov, muv, piv, vp, vs, y, t, de, zplat1, zplat2, rplat;
	float *pts, ts, tp, sumu, sumpi, ws, *ri, *d, *ti, *dl, *vsl, z, r;
	float **checkp, **checks, **checkrho; 
	int   i, j, l, ii, jj, nk, k, nl;
	char filename_mu[STRING_SIZE];
	char filename_rho[STRING_SIZE]; 
	char filename_pi[STRING_SIZE];
				
	sprintf(filename_mu,"%s.mu",MFILE);
	sprintf(filename_rho,"%s.rho",MFILE);
	sprintf(filename_pi,"%s.pi",MFILE);
	
	
	/*-----------------------------------------------------------------------*/

	/* loop over global grid */
	for (i=1;i<=NXG;i++){
            for (j=1;j<=NYG;j++){
	                
	                  y = j*DH;
	                
			  vs = 1150.0;
			  vp = 2000.0;
			  rhov = 1800.0;
			  
			  if(y<0.01){
			   
			    vs = 2310.0;  
			    vp = 4000.0;  
			  rhov = 1800.0;
			   
			  }
	    
	                   
	    
			/* only the PE which belongs to the current global gridpoint 
			is saving model parameters in his local arrays */
			if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
				ii = i-POS[1]*NX;
				jj = j-POS[2]*NY;
				
				u[jj][ii]    = vs;
				rho[jj][ii]  = rhov;
				pi[jj][ii]   = vp;
			}
			

		}
	}	

	/* each PE writes his model to disk */
	writemod(filename_rho,rho,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(filename_rho,3);
	
	/* each PE writes his model to disk */
	writemod(filename_pi,pi,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(filename_pi,3);
	                        
	/* each PE writes his model to disk */
	writemod(filename_mu,u,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(filename_mu,3);
	                        
}


