/*
 *   Elastic homogeneous, anisotropic TTI zinc model
 * 
 *   from 
 *
 *   S. Operto, J. Virieux, A. Ribodetti, and J. E. Anderson. 
 *   Finite-difference frequency-domain modeling of visco-acoustic 
 *   wave propagation in two-dimensional TTI media. 
 *   Geophysics, 74 (5):T75--T95, 2009
 *   
 *   last update 05/02/2017 
 *   D. Koehn
 *
 */

#include "fd.h"

void model_elastic_TTI(float  **  rho, float **  c11, float **  c13, float **  c33, float **  c44, float **  theta){

	/*--------------------------------------------------------------------------*/
	/* global variables */
	extern int NX, NY, NXG, NYG,  POS[3], MYID;
	extern char  MFILE[STRING_SIZE], INV_MODELFILE[STRING_SIZE];	
	extern float DH;
	
        /* local variables */	
	int i, j, ii, jj;
	char filename[STRING_SIZE];
	
	/* anisotropic parameters */
        float vp0 = 2955.06, vsv = 2361.67, rhoh=7100.0, delta = 2.70968, epsilon = 0.830645, thetah = 45.0;

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
					theta[jj][ii]=thetah * M_PI / 180.0;
					rho[jj][ii]=rhoh;
					
				}
			}
		}	

	/* each PE writes his model to disk and PE 0 merges model files */
	sprintf(filename,"%s.denise.c11",MFILE);
	writemod(filename,c11,3);
	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(filename,3);
	
	sprintf(filename,"%s.denise.c13",MFILE);
        writemod(filename,c13,3);
	MPI_Barrier(MPI_COMM_WORLD);
	                           
	if (MYID==0) mergemod(filename,3);
	
	sprintf(filename,"%s.denise.c33",MFILE);
	writemod(filename,c33,3);
	MPI_Barrier(MPI_COMM_WORLD);
	                        
	if (MYID==0) mergemod(filename,3);

        sprintf(filename,"%s.denise.c44",MFILE);
	writemod(filename,c44,3);
	MPI_Barrier(MPI_COMM_WORLD);
	                        
	if (MYID==0) mergemod(filename,3);

	sprintf(filename,"%s.denise.theta",MFILE);
	writemod(filename,theta,3);
	MPI_Barrier(MPI_COMM_WORLD);
	                        
	if (MYID==0) mergemod(filename,3);

        /* clean up temporary files */
        MPI_Barrier(MPI_COMM_WORLD);

        sprintf(filename,"%s.denise.c11.%i%i",MFILE,POS[1],POS[2]);
        remove(filename);

        sprintf(filename,"%s.denise.c13.%i%i",MFILE,POS[1],POS[2]);
        remove(filename);

        sprintf(filename,"%s.denise.c33.%i%i",MFILE,POS[1],POS[2]);
        remove(filename);

        sprintf(filename,"%s.denise.c44.%i%i",MFILE,POS[1],POS[2]);
        remove(filename);

        sprintf(filename,"%s.denise.theta.%i%i",MFILE,POS[1],POS[2]);
        remove(filename);

}



