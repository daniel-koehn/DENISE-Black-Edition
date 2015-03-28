/*
 *   homogeneous half space 
 *   last update 05.02.99, T. Bohlen
 */

#include "fd.h"

double drand48(void);
void srand48(long seedval);

void model(float  **  rho, float **  pi, float **  u, 
float **  taus, float **  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, *FL, TAU;
	extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	
	/* local variables */
	float rhov, muv, piv, vp, vs, r;
	float *pts, ts, tp, sumu, sumpi, ws;
	int i, j, l, ii, jj;
//	FILE *fp;
	char modfile[STRING_SIZE];
	

	/* paramemeters for computation of density and Vs */
	const float vp0=2800.0, vpfluc=0.0;

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


	srand48(NXG*NYG);

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
				r=(float)drand48();
				vp=vp0+(2.0*r-1.0)*vpfluc;
				
				vs=-314.59 + 0.61*vp;   
				rhov=1498.0 + 0.22*vp;

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
					u[jj][ii]=muv;
					rho[jj][ii]=1.0/rhov;
					pi[jj][ii]=piv;
				}
			}
		}	

	sprintf(modfile,"model/model.bin");

	/* each PE writes his model to disk */
	writemod(modfile,pi,3);

	MPI_Barrier(MPI_COMM_WORLD);

	/*  model files are then merged into one file by PE 0 */
	if (MYID==0) mergemod(modfile,3);

	free_vector(pts,1,L);
}



