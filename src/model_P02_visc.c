/*
 *   Model homogeneous half space
 *   last update 11.04.02, T. Bohlen
 */

#include "fd.h"

void model(float  **  rho, float **  pi, float **  u, float **  taus, float **  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	extern char  MFILE[STRING_SIZE];	
	extern char INV_MODELFILE[STRING_SIZE];
	extern float DH, *FL, TAU, DT;
		/* local variables */
	float vp, vs, rhov, grad1, grad2, grad3, y, ts, tp, muv, piv, *pts, lkappa, lmu;
	float Qp, Qs, Qpinv, llambda;
	int i, j, ii, jj, l;
	char modfile[STRING_SIZE]; 
	
	/* parameters for homogenous half-space */
	/* perfect model parameters */
	const float vp1=2715.0, vs1=1677.0, rho1=2120.0, Qmu=20.0, Qkappa=10000.0;
	/*const float vp1=2700.0, vs1=1970.0, rho1=1190.0, Qmu=50.0, Qkappa=10000.0;*/
	
	/*-----------------------------------------------------------------------*/
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
	        eta[l]=DT/pts[l];
	}	
	
	
	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
						
				vp=vp1;
				vs=vs1;
				rhov=rho1;
				
				/* calculate kappa and mu*/
				/*lmu=vs*vs*rhov;
				llambda = rhov*(vp*vp*rhov - 2.0*lmu); 
				lkappa = llambda + ((2.0/3.0)*lmu);
			        Qpinv = (1.0/Qmu) + (lkappa/(lkappa+(4.0*lmu/3.0)))*((1.0/Qkappa)-(1.0/Qmu));
				Qp = 1.0/Qpinv;*/
				
				Qp = Qkappa;
				Qs = Qmu;
								
				ts=2.0/Qs;
				tp=2.0/Qp;
				
				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;

					u[jj][ii]=vs;
					rho[jj][ii]=rhov;
					pi[jj][ii]=vp;
					taus[jj][ii]=ts;
					taup[jj][ii]=tp;
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

free_vector(pts,1,L);
}



