/*------------------------------------------------------------------------
 *  Read VTI model properties (c11, c13, c33, c44, density) from external files  
 *
 *  D. Koehn
 *  Kiel, 02.02.2017
 *  ----------------------------------------------------------------------*/


#include "fd.h"

void readmod_elastic_VTI(float  **  rho, float **  c11, float **  c13, float **  c33, float **  c44){

	extern int NX, NY, NXG, NYG,  POS[3], MYID, INVMAT1;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

		
	/* local variables */
	float rhov, c11v, c13v, c33v, c44v;
	int i, j, ii, jj;
	FILE *fp_c11, *fp_c13, *fp_c33, *fp_c44, *fp_rho;
	char filename[STRING_SIZE];


	fprintf(FP,"\n...reading model information from external model-files...\n");
           
	   
	/* read density and VTI elastic tensor parameters */
        /* ---------------------------------------------- */

	fprintf(FP,"\t c11:\n\t %s.c11\n\n",MFILE);
	sprintf(filename,"%s.c11",MFILE);
	fp_c11=fopen(filename,"r");
	if (fp_c11==NULL) err(" Could not open model file for c11 ! ");

	fprintf(FP,"\t c13:\n\t %s.c13\n\n",MFILE);
	sprintf(filename,"%s.c13",MFILE);
	fp_c13=fopen(filename,"r");
	if (fp_c13==NULL) err(" Could not open model file for c13 ! ");

	fprintf(FP,"\t c33:\n\t %s.c33\n\n",MFILE);
	sprintf(filename,"%s.c33",MFILE);
	fp_c33=fopen(filename,"r");
	if (fp_c33==NULL) err(" Could not open model file for c33 ! ");

	fprintf(FP,"\t c44:\n\t %s.c44\n\n",MFILE);
	sprintf(filename,"%s.c44",MFILE);
	fp_c44=fopen(filename,"r");
	if (fp_c44==NULL) err(" Could not open model file for c44 ! ");

	fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
	sprintf(filename,"%s.rho",MFILE);
	fp_rho=fopen(filename,"r");
	if (fp_rho==NULL) err(" Could not open model file for densities ! ");
	   

	/* loop over global grid */
	for (i=1;i<=NXG;i++){
	    for (j=1;j<=NYG;j++){

	        fread(&c11v, sizeof(float), 1, fp_c11);
	        fread(&c13v, sizeof(float), 1, fp_c13);
	        fread(&c33v, sizeof(float), 1, fp_c33);
	        fread(&c44v, sizeof(float), 1, fp_c44);
		fread(&rhov, sizeof(float), 1, fp_rho);
				
		/* only the PE which belongs to the current global gridpoint 
		is saving model parameters in his local arrays */
		if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){

		    ii=i-POS[1]*NX;
		    jj=j-POS[2]*NY;
    		    
		    c11[jj][ii]=c11v;
		    c13[jj][ii]=c13v;
		    c33[jj][ii]=c33v;
		    c44[jj][ii]=c44v;
                    rho[jj][ii]=rhov;                                                                    
				
		}
	    }
	}
	
	fclose(fp_c11);
	fclose(fp_c13);
	fclose(fp_c33);
	fclose(fp_c44);
	fclose(fp_rho);
	
	
	/* each PE writes his model to disk */
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

}




