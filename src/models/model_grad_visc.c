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
	float vp, vs, rhov, grad1, grad2, grad3, y, ts, tp, muv, piv, *pts;
	int i, j, ii, jj, l;
	char modfile[STRING_SIZE]; 
	
	/* parameters for layer 1 velocity in m/s, depth h in meter */
	const float vp1=600.0, vs1=280.0, rho1=1800.0, h=10.0;
	
	/* parameters for layer 2 due to calculation of grad1, grad2 and grad3*/
	const float vp2=2500.0, vs2=800.0, rho2=2000.0;
	
	/*-----------------------------------------------------------------------*/
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
	        eta[l]=DT/pts[l];
	}
	
	
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
				
				ts=TAU;
				tp=TAU;
				
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



