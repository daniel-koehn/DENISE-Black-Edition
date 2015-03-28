/*------------------------------------------------------------------------
 *   calculate test step length for material parameter update
 *   
 *   Daniel Koehn
 *   last update 9.11.2007
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
void model_freq_out(float ** ppi, float  **  rho, float **  pu, int iter, float freq){


	/*--------------------------------------------------------------------------*/
	FILE *FP1;
	/* extern variables */
	extern int POS[3], MYID;
	extern char INV_MODELFILE[STRING_SIZE];
	

	/* local variables */
        char modfile[STRING_SIZE];

        sprintf(modfile,"%s_vp_stage_%d.bin",INV_MODELFILE,iter);
	writemod(modfile,ppi,3);
	MPI_Barrier(MPI_COMM_WORLD);
                                                                                                        
	if (MYID==0) mergemod(modfile,3);
	MPI_Barrier(MPI_COMM_WORLD); 
	sprintf(modfile,"%s_vp_stage_%d.bin.%i%i",INV_MODELFILE,iter,POS[1],POS[2]);
	remove(modfile);                                                                                                                        
                          
	sprintf(modfile,"%s_vs_stage_%d.bin",INV_MODELFILE,iter);
	writemod(modfile,pu,3);
	MPI_Barrier(MPI_COMM_WORLD);
                                                                                                        
	if (MYID==0) mergemod(modfile,3);
	MPI_Barrier(MPI_COMM_WORLD); 
	sprintf(modfile,"%s_vs_stage_%d.bin.%i%i",INV_MODELFILE,iter,POS[1],POS[2]);
	remove(modfile);                                                                                                                        
                                                                                                                                
	sprintf(modfile,"%s_rho_stage_%d.bin",INV_MODELFILE,iter);
	writemod(modfile,rho,3);
	MPI_Barrier(MPI_COMM_WORLD);
                                                                                                                                                                        
	if (MYID==0) mergemod(modfile,3);
	MPI_Barrier(MPI_COMM_WORLD); 
	sprintf(modfile,"%s_rho_stage_%d.bin.%i%i",INV_MODELFILE,iter,POS[1],POS[2]);
	remove(modfile);

}

