/*------------------------------------------------------------------------
 *  Read viscoelastic model properties (vs,density,Qs) from files  
 *
 *  D. Koehn
 *  Kiel, 03.12.2017
 *  ----------------------------------------------------------------------*/


#include "fd.h"

void readmod_visc_SH(float  **  rho, float **  u, float ** taus, float * eta){

	extern float DT, *FL;
	extern int L;
	extern int NX, NY, NXG, NYG,  POS[3], MYID, INVMAT1;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

		
	/* local variables */
	float rhov, muv, vs, qs, *pts;
	int i, j, ii, jj, l;
	FILE *fp_vs, *fp_rho, *fp_qs;
	char filename[STRING_SIZE];
	
	
	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}


	   fprintf(FP,"\n...reading model information from model-files...\n");
           
	   /* read density and seismic velocities */
	   /* ----------------------------------- */
	   if(INVMAT1==1){ 

	   fprintf(FP,"\t Vs:\n\t %s.vs\n\n",MFILE);
	   sprintf(filename,"%s.vs",MFILE);
	   fp_vs=fopen(filename,"r");
	   if (fp_vs==NULL) err(" Could not open model file for Vs ! ");

	   fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
	   sprintf(filename,"%s.rho",MFILE);
	   fp_rho=fopen(filename,"r");
	   if (fp_rho==NULL) err(" Could not open model file for densities ! ");	   
	   
	   fprintf(FP,"\t Qs:\n\t %s.qs\n\n",MFILE);
	   sprintf(filename,"%s.qs",MFILE);
	   fp_qs=fopen(filename,"r");
	   if (fp_qs==NULL) err(" Could not open model file for Qs-values ! ");
           }
	   
	   /* read density and Lame parameters */
	   /* ----------------------------------- */
	   if(INVMAT1==3){ 

	   fprintf(FP,"\t Lame parameter mu:\n\t %s.mu\n\n",MFILE);
	   sprintf(filename,"%s.mu",MFILE);
	   fp_vs=fopen(filename,"r");
	   if (fp_vs==NULL) err(" Could not open model file for Lame parameter mu ! ");

	   fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
	   sprintf(filename,"%s.rho",MFILE);
	   fp_rho=fopen(filename,"r");
	   if (fp_rho==NULL) err(" Could not open model file for densities ! ");	   
	   
	   fprintf(FP,"\t Qs:\n\t %s.qs\n\n",MFILE);
	   sprintf(filename,"%s.qs",MFILE);
	   fp_qs=fopen(filename,"r");
	   if (fp_qs==NULL) err(" Could not open model file for Qs-values ! ");
           }
	   

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			fread(&vs, sizeof(float), 1, fp_vs);
			fread(&rhov, sizeof(float), 1, fp_rho);
			fread(&qs, sizeof(float), 1, fp_qs);
				
			/* only the PE which belongs to the current global gridpoint 
			is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
                                
				u[jj][ii]=vs;
                                rho[jj][ii]=rhov;
				taus[jj][ii]=2.0/qs;
				
				}
			}
		}
	

	fclose(fp_vs);
	fclose(fp_rho);
	fclose(fp_qs);
	
	
	/* each PE writes his model to disk */
	
	sprintf(filename,"%s.denise.mu",MFILE);
        writemod(filename,u,3);
	MPI_Barrier(MPI_COMM_WORLD);
	                           
	if (MYID==0) mergemod(filename,3);
	
	sprintf(filename,"%s.denise.rho",MFILE);
	writemod(filename,rho,3);
	MPI_Barrier(MPI_COMM_WORLD);
	                        
	if (MYID==0) mergemod(filename,3);

	free_vector(pts,1,L);
}




