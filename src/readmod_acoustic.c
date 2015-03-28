/* $Id: readmod_acoustic.c,v 1.1.1.1 2007/12/28 18:25:05 koehn Exp $*/
/*------------------------------------------------------------------------
 *   Read elastic model properties (vp,vs,density) from files  
 *
 *  ----------------------------------------------------------------------*/


/* This file contains function readmod, which has the purpose
   to read data from model-files for viscoelastic simulation */

#include "fd.h"

void readmod_acoustic(float  **  rho, float **  pi, float **  u){

	extern int NX, NY, NXG, NYG,  POS[3], MYID;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;
	extern float RHOHOM;
		
	/* local variables */
	float rhov, muv, piv, vp, vs;
	int i, j, ii, jj;
	FILE *fp_vs, *fp_vp, *fp_rho;
	char filename[STRING_SIZE];




	   fprintf(FP,"\n...reading modell information from modell-files...\n");

	   fprintf(FP,"\t P-wave velocities:\n\t %s\n\n",MFILE);
	   sprintf(filename,"%s",MFILE);
	   fp_vp=fopen(filename,"r");
	   if (fp_vp==NULL) err(" Could not open modell file for P velocities ! ");

/*	   fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
	   sprintf(filename,"%s.rho",MFILE);
	   fp_rho=fopen(filename,"r");
	   if (fp_rho==NULL) err(" Could not open modell file for densities ! ");
*/
	   

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			fread(&vp, sizeof(float), 1, fp_vp);
			/*
			fread(&vs, sizeof(float), 1, fp_vs);
			fread(&rhov, sizeof(float), 1, fp_rho);
			*/
			/* Gardner */
			vs = 1e-5;
			rhov = RHOHOM;

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
	




	fclose(fp_vp);
	/*
	fclose(fp_vs);
	fclose(fp_rho);
	*/
	
	
	/* each PE writes his model to disk */
	   
	   
	sprintf(filename,"%s.pi",MFILE);
	writemod(filename,pi,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(filename,3);

	sprintf(filename,"%s.rho",MFILE);
	writemod(filename,rho,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(filename,3);

	sprintf(filename,"%s.u",MFILE);
	writemod(filename,u,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(filename,3);

}




