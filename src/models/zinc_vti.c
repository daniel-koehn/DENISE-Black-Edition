/*
 *   Elastic homogeneous, anisotropic VTI zinc model
 * 
 *   from 
 *
 *   S. Operto, J. Virieux, A. Ribodetti, and J. E. Anderson. 
 *   Finite-difference frequency-domain modeling of visco-acoustic 
 *   wave propagation in two-dimensional TTI media. 
 *   Geophysics, 74 (5):T75--T95, 2009
 *   
 *   last update 03/02/2017 
 *   D. Koehn
 *
 */

#include "fd.h"

void model_elastic_VTI(float  **  rho, float **  c11, float **  c13, float **  c33, float **  c44){

	/*--------------------------------------------------------------------------*/
	/* global variables */
	extern int NX, NY, NXG, NYG,  POS[3], MYID;
	extern char  MFILE[STRING_SIZE], INV_MODELFILE[STRING_SIZE];	
	extern float DH;
	
        /* local variables */	
	int i, j, ii, jj;
	char modfile[STRING_SIZE];
	
	/* anisotropic parameters */
        float vp0 = 2955.06, vsv = 2361.67, rhoh=7100.0, delta = 2.70968, epsilon = 0.830645;

        /* transform Thomsen's parameters to elastic tensor components */
        float c33h = rhoh * vp0 * vp0;
        float c44h = rhoh * vsv * vsv;
	float c11h = c33h * (1 + 2.0 * epsilon);
	float c13h = sqrt((c33h-c44h) * (c33h-c44h) + 2.0 * delta * c33h * (c33h - c44h)) - c44h;
		
	/*-----------------------------------------------------------------------*/

		
	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
											
       				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;

					c11[jj][ii]=c11h;
					c13[jj][ii]=c13h;
					c33[jj][ii]=c33h;
					c44[jj][ii]=c44h;
					rho[jj][ii]=rhoh;
					
				}
			}
		}	

/*sprintf(modfile,"%s_rho_it_0.bin",INV_MODELFILE);
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
if (MYID==0) mergemod(modfile,3);*/		

}



