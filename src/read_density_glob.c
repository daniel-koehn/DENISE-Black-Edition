/*------------------------------------------------------------------------
 *   Read global density model from files  
 *
 *  D. Koehn
 *  Kiel, the 7th of November 2014
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void read_density_glob(float  ** rho_grav, int sws){

	extern int NXG, NYG, MYID, INVMAT1;
	extern char  MFILE[STRING_SIZE];
	extern char  JACOBIAN[STRING_SIZE];	

	/* local variables */
	float rhov;
	int i, j, ii, jj;
	FILE *fp_rho;
	char filename[STRING_SIZE];

        /* read initial model file */
        if(sws==1){
	  sprintf(filename,"%s.rho",MFILE);
	}
	
	/* read density model from current iteration */
	if(sws==2){
	  sprintf(filename,"%s_tmp.rho",JACOBIAN);
	  /*printf("%s",filename);*/
	}
	
	fp_rho=fopen(filename,"r");
	if (fp_rho==NULL) err(" Could not open model file for densities ! ");
      
	/* loop over global grid */
	for (i=1;i<=NXG;i++){
		for (j=1;j<=NYG;j++){
	
			fread(&rhov, sizeof(float), 1, fp_rho);
			rho_grav[j][i]=rhov;
			
		}
	}
	
	fclose(fp_rho);
	
}




