/*
 * Apply median filter at source positions
 *
 * Daniel Koehn
 * Kiel, the 31st of July 2013
 */

#include "fd.h"

void median_src(float ** waveconv,float ** taper_coeff, float **srcpos, int nshots, int **recpos, int ntr, int iter, int sws)
{

	/* extern variables */
        extern float DH;
	extern int FREE_SURF, NX, NY, NXG, NYG, IDX, IDY;
	extern int NPROCX, NPROCY, MYID, POS[3];
	extern char JACOBIAN[STRING_SIZE];
	
	/* local variables */
	int i, j, h, ii, jj, n;
        char jac[STRING_SIZE];
        float gradtmp;
        float ** gradtmp1;
        FILE *FP3, *fp_grad;        

        /* define temporary gradient matrix */
        gradtmp1 = matrix(1,NYG,1,NXG);        

        /* temporarily save gradient for median filtering at source positions */
        sprintf(jac,"%s_wavenumber.old.%i%i",JACOBIAN,POS[1],POS[2]);
        FP3=fopen(jac,"wb");
                        
        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                fwrite(&waveconv[j][i],sizeof(float),1,FP3);
            }                       
        }                           
                                    
        fclose(FP3);                
                                    
        MPI_Barrier(MPI_COMM_WORLD);
                                  
        /* merge gradient file */ 
        sprintf(jac,"%s_wavenumber.old",JACOBIAN);
        if (MYID==0) mergemod(jac,3);

        MPI_Barrier(MPI_COMM_WORLD);
        
	fp_grad=fopen(jac,"rb");
	
	if (fp_grad==NULL) err(" Could not open gradient file ! ");
	
	/* load merged gradient */
	for (i=1;i<=NXG;i++){
	   for (j=1;j<=NYG;j++){
	        
	            fread(&gradtmp, sizeof(float), 1, fp_grad);
	            gradtmp1[j][i] = gradtmp;
			
            }
	}
	
	fclose(fp_grad);        

	/* apply median filter at shot positions */
        for (n=1;n<=nshots;n++){
          i = iround(srcpos[1][n]/DH);
          j = iround(srcpos[2][n]/DH);
          gradtmp1[j][i] = (gradtmp1[j-1][i] + gradtmp1[j+1][i] + gradtmp1[j][i-1] + gradtmp1[j][i+1] + gradtmp1[j-1][i-1] + gradtmp1[j-1][i+1] + gradtmp1[j+1][i+1] + gradtmp1[j+1][i-1])/8.0;        
        }

         /* distribute spatial filtered gradient on computational nodes */
	 for (i=1;i<=NXG;i++){
	    for (j=1;j<=NYG;j++){

			if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				waveconv[jj][ii]=gradtmp1[j][i];
			}
			
	     }
	  }

          /* clean up memory */
          free_matrix(gradtmp1,1,NYG,1,NXG);


}



