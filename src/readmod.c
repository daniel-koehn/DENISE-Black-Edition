/*------------------------------------------------------------------------
 *   Read viscoelastic model properties (vp,vs,density,Qp,Qs) from files  
 *
 *  Copyright (c)  T. Bohlen
 *  ----------------------------------------------------------------------*/


/* This file contains function readmod, which has the purpose
   to read data from model-files for viscoelastic simulation */

#include "fd.h"

void readmod(float  **  rho, float **  pi, float **  u, float ** taus, float ** taup, float * eta){

	extern float DT, *FL;
	extern int L;
	extern int NX, NY, NXG, NYG,  POS[3], MYID, INVMAT1;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

		
	/* local variables */
	float rhov, muv, piv, vp, vs, qp, qs, *pts;
	int i, j, ii, jj, l;
	FILE *fp_vs, *fp_vp, *fp_rho, *fp_qp, *fp_qs;
	char filename[STRING_SIZE];
	
	
	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}





	   fprintf(FP,"\n...reading model information from modell-files...\n");
           
	   /* read density and seismic velocities */
	   /* ----------------------------------- */
	   if(INVMAT1==1){ 
	   fprintf(FP,"\t Vp:\n\t %s.vp\n\n",MFILE);
	   sprintf(filename,"%s.vp",MFILE);
	   fp_vp=fopen(filename,"r");
	   if (fp_vp==NULL) err(" Could not open model file for Vp ! ");


	   fprintf(FP,"\t Vs:\n\t %s.vs\n\n",MFILE);
	   sprintf(filename,"%s.vs",MFILE);
	   fp_vs=fopen(filename,"r");
	   if (fp_vs==NULL) err(" Could not open model file for Vs ! ");

	   fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
	   sprintf(filename,"%s.rho",MFILE);
	   fp_rho=fopen(filename,"r");
	   if (fp_rho==NULL) err(" Could not open model file for densities ! ");
	   
	   fprintf(FP,"\t Qp:\n\t %s.qp\n\n",MFILE);
	   sprintf(filename,"%s.qp",MFILE);
	   fp_qp=fopen(filename,"r");
	   if (fp_qp==NULL) err(" Could not open model file for Qp-values ! ");
	   
	   fprintf(FP,"\t Qs:\n\t %s.qs\n\n",MFILE);
	   sprintf(filename,"%s.qs",MFILE);
	   fp_qs=fopen(filename,"r");
	   if (fp_qs==NULL) err(" Could not open model file for Qs-values ! ");
           }
	   
	   /* read density and Lame parameters */
	   /* ----------------------------------- */
	   if(INVMAT1==3){ 
	   fprintf(FP,"\t Lame parameter lambda:\n\t %s.lam\n\n",MFILE);
	   sprintf(filename,"%s.lam",MFILE);
	   fp_vp=fopen(filename,"r");
	   if (fp_vp==NULL) err(" Could not open model file for Lame parameter lambda ! ");


	   fprintf(FP,"\t Lame parameter mu:\n\t %s.mu\n\n",MFILE);
	   sprintf(filename,"%s.mu",MFILE);
	   fp_vs=fopen(filename,"r");
	   if (fp_vs==NULL) err(" Could not open model file for Lame parameter mu ! ");

	   fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
	   sprintf(filename,"%s.rho",MFILE);
	   fp_rho=fopen(filename,"r");
	   if (fp_rho==NULL) err(" Could not open model file for densities ! ");
	   
	   fprintf(FP,"\t Qp:\n\t %s.qp\n\n",MFILE);
	   sprintf(filename,"%s.qp",MFILE);
	   fp_qp=fopen(filename,"r");
	   if (fp_qp==NULL) err(" Could not open model file for Qp-values ! ");
	   
	   fprintf(FP,"\t Qs:\n\t %s.qs\n\n",MFILE);
	   sprintf(filename,"%s.qs",MFILE);
	   fp_qs=fopen(filename,"r");
	   if (fp_qs==NULL) err(" Could not open model file for Qs-values ! ");
           }
	   

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			fread(&vp, sizeof(float), 1, fp_vp);
			fread(&vs, sizeof(float), 1, fp_vs);
			fread(&rhov, sizeof(float), 1, fp_rho);
			fread(&qp, sizeof(float), 1, fp_qp);
			fread(&qs, sizeof(float), 1, fp_qs);
			
				
			/* only the PE which belongs to the current global gridpoint 
			is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
                                
				u[jj][ii]=vs;
                                rho[jj][ii]=rhov;
                                pi[jj][ii]=vp;
				taus[jj][ii]=2.0/qs;
				taup[jj][ii]=2.0/qp;
				
				}
			}
		}
	




	fclose(fp_vp);
	fclose(fp_vs);
	fclose(fp_rho);
	fclose(fp_qp);
	fclose(fp_qs);
	
	
	/* each PE writes his model to disk */
	sprintf(filename,"%s.fdveps.pi",MFILE);
	writemod(filename,pi,3);
	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(filename,3);
	
	sprintf(filename,"%s.fdveps.mu",MFILE);
        writemod(filename,u,3);
	MPI_Barrier(MPI_COMM_WORLD);
	                           
	if (MYID==0) mergemod(filename,3);
	
	sprintf(filename,"%s.fdveps.rho",MFILE);
	writemod(filename,rho,3);
	MPI_Barrier(MPI_COMM_WORLD);
	                        
	if (MYID==0) mergemod(filename,3);

	free_vector(pts,1,L);
}




