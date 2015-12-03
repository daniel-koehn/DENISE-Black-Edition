/*------------------------------------------------------------------------
 *   Read background density for gravity modelling from file  
 *
 *  Copyright (c)  D. Wehner
 *  last update 13.10.2015
 *  ----------------------------------------------------------------------*/


#include "fd.h"

void read_back_density(float  **  rho_back){

	extern int NXG, NYG, MYID, BACK_DENSITY, NX, NY, POS[3];
	extern char  MFILE[STRING_SIZE], DFILE[STRING_SIZE];
	extern FILE *FP;	

	/* local variables */
	float rhov;
	int i, j, ii, jj;
	FILE *fp_rho;
	char filename[STRING_SIZE];

	
	/* read background density from initial model*/
	/* ----------------------------------- */	
	if(BACK_DENSITY == 1){
		fprintf(FP,"\n...reading background density information from initial model-file...\n");
	
		sprintf(filename,"%s.rho",MFILE);
		fp_rho=fopen(filename,"r");
		if (fp_rho==NULL) err(" Could not open model file for density ! ");
		
		/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
	
				fread(&rhov, sizeof(float), 1, fp_rho);					
				rho_back[j][i]=rhov;
			
			}
		}
	
	fclose(fp_rho);
		
		
	}	
	

	/* read background density from initial model*/
	/* ----------------------------------- */
	/*sprintf(filename,"%s.rho",DFILE);*/
	if(BACK_DENSITY == 2){
		fprintf(FP,"\n...reading background density information from self-defined model-file...");
		
		fp_rho=fopen(DFILE,"r");
		if (fp_rho==NULL) {err(" Could not open model file for density ! ");}
        
		/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
	
				fread(&rhov, sizeof(float), 1, fp_rho);
				rho_back[j][i]=rhov;				
			
			}
		}
	
	fclose(fp_rho);
	
	}
		
	fprintf(FP,"...finished...\n");

}
