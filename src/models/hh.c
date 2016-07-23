/*
 *   Model layer over half space (homogeneous)
 *   update 11.04.02, T. Bohlen
 *   last update 11. Okt. 2011 M. Schaefer
 */

#include "fd.h"

void model(float  **  rho, float **  pi, float **  u, 
float **  taus, float **  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, *FL, TAU, DH;
	extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	extern char  MFILE[STRING_SIZE];	

		/* local variables */
	float rhov, muv, piv, vp, vs, y;
	float *pts, ts, tp, sumu, sumpi, ws;
	int i, j, l, ii, jj;
	 
	
	/* parameters for layer 1 */
	const float vp1=450.0, vs1=270.0, rho1=1800.0, h=5.0;
	
	/* parameters for layer 2 */
	const float vp2=450.0, vs2=270.0, rho2=1800.0;
	
	
	/*-----------------------------------------------------------------------*/


	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}

	ts=TAU;  
	tp=TAU;

	ws=2.0*PI*FL[1];
	
	sumu=0.0; 
	sumpi=0.0;
	for (l=1;l<=L;l++){
		sumu=sumu+((ws*ws*pts[l]*pts[l]*ts)/(1.0+ws*ws*pts[l]*pts[l]));
		sumpi=sumpi+((ws*ws*pts[l]*pts[l]*tp)/(1.0+ws*ws*pts[l]*pts[l]));
	}

		

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			
				y=(float)j*DH;
			
				if (y<=h){
				 vp=vp1; vs=vs1; rhov=rho1; }

				
				 else{
 				 vp=vp2; vs=vs2; rhov=rho2; }
                    
				
				muv=vs*vs*rhov/(1.0+sumu);
				piv=vp*vp*rhov/(1.0+sumpi);

				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;

					taus[jj][ii]=ts;
					taup[jj][ii]=tp;
					u[jj][ii]=vs1;
					rho[jj][ii]=rho1;
					pi[jj][ii]=vp1;
				}
			}
		}	

		

	
	/* each PE writes his model to disk */

	writemod(MFILE,pi,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(MFILE,3);

	free_vector(pts,1,L);
}



