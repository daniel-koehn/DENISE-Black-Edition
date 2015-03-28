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
	extern float DH;	

		/* local variables */
	float muv, piv, Vp, Vs, Rho, y, x, z ;
	int i, j, ii, jj;
	char fname[STRING_SIZE];

	
	/* parameters for layer 1 */
	const float vp1=320.0, vs1=200.0, rho1=1700.0, z1=2.0;
	const float vp2=320.0, vs2=200.0, rho2=1800.0, z2=20.0;
	const float vp3=1550.0, vs3=250.0, rho3=2000.0;
	const float w=2.0*PI/100.0, a=3.0;
	
	
	/*-----------------------------------------------------------------------*/



		

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			
			y=(float)j*DH;
			x=(float)i*DH;
			z=z2+a*sin(w*x);
			if (y<=z1){Vp=vp1;Vs=vs1;Rho=rho1;}
			else if (y<=z){Vp=vp2;Vs=vs2;Rho=rho2;}
			else {Vp=vp3;Vs=vs3;Rho=rho3;}
			

				
       				muv=Vs*Vs*Rho;
				piv=Vp*Vp*Rho;

				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;

					u[jj][ii]=muv;
					rho[jj][ii]=Rho;
					pi[jj][ii]=piv;
				}
			}
		}	

		
	/* each PE writes his model to disk */
	sprintf(fname,"%s.mu",MFILE);	        
	writemod(fname,u,3);

	MPI_Barrier(MPI_COMM_WORLD);

	/*  model files are then merged into one file by PE 0 */
	if (MYID==0) mergemod(fname,3);
}



