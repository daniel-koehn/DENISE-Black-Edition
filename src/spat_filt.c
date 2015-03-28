/*
 * taper gradient with a gaussian frame to damp inversion artefacts near the sources and receivers
 sws == 1 vertical taper (for tomography geometry)
 sws == 2 horizontal taper (for reflection geometry)
 sws == 3 local radial taper at the source and receiver positions
 sws == 4  
 */

#include "fd.h"

void spat_filt(float ** waveconv, int iter, int sws)
{

	/* extern variables */

        extern float DH;
	extern int FREE_SURF, NX, NY, NXG, NYG;
	extern int NPROCX, NPROCY, MYID, POS[3];
	extern char JACOBIAN[STRING_SIZE];
	extern FILE *FP;
	
	/* local variables */
	int i, j, h, ifw, ii, jj, n, xb, yb, xe, ye, taperlength,taperlength2;
	int ijc, iy, ix, iii, jjj, xx, yy, srctaper_gridpt, i1, j1, filtsize;

	float **waveconvtmp, **waveconvtmps, grad, normgauss;
	char jac[STRING_SIZE];
	FILE *fp_grad;
	
	extern int SPAT_FILT_SIZE, SPAT_FILT_1, SPAT_FILT_ITER;
	/*SPAT_FILT_SIZE=23;
	SPAT_FILT_1=45;*/
	
	if (SPAT_FILT_ITER){   /* filter length is reduced in each iteration */
		SPAT_FILT_SIZE=SPAT_FILT_SIZE-iter+1;}
	
	/*if(SPAT_FILT_SIZE>=iter){*/ /* application condition for the spatial filter */
	if(SPAT_FILT_SIZE>=1){
	
        waveconvtmp = matrix(1,NYG,1,NXG);
	waveconvtmps = matrix(1,NYG,1,NXG);
	
	if(sws==1){
	sprintf(jac,"%s_g.old",JACOBIAN);}   /*sprintf(jac,"%s_g.it%i",JACOBIAN,iter);}*/

	if(sws==2){
	sprintf(jac,"%s_g_u.old",JACOBIAN);}    /*sprintf(jac,"%s_g_u.it%i",JACOBIAN,iter);}*/
	
	if(sws==3){
	sprintf(jac,"%s_g_rho.old",JACOBIAN);}    /*sprintf(jac,"%s_g_rho.it%i",JACOBIAN,iter);}*/
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(MYID==0){
	fp_grad=fopen(jac,"rb");
	if (fp_grad==NULL) err(" Could not open gradient file ! ");
	
	/* load merged gradient */
	for (i=1;i<=NXG;i++){
	   for (j=1;j<=NYG;j++){
                        fread(&grad, sizeof(float), 1, fp_grad);
			waveconvtmp[j][i]=grad;
			waveconvtmps[j][i]=grad;
		}
	}
	
	fclose(fp_grad);
	
	/* apply spatial Gaussian filter on gradient */
	for (i=1+SPAT_FILT_SIZE;i<=NXG-SPAT_FILT_SIZE;i++){
	   for (j=SPAT_FILT_1;j<=NYG-SPAT_FILT_SIZE;j++){

                     waveconvtmps[j][i]=0.0;   
		     normgauss=0.0;
		     
		     for (i1=i-SPAT_FILT_SIZE;i1<=i+SPAT_FILT_SIZE;i1++){
	                 /*for (j1=j+SPAT_FILT_SIZE;j1<=j+SPAT_FILT_SIZE;j1++){*/
			 
                         waveconvtmps[j][i] += waveconvtmp[j][i1] * (1/(sqrt(2*PI)*SPAT_FILT_SIZE)) * exp((-1/2)*((sqrt(((i-i1)*(i-i1))+((i-i1)*(i-i1))))*(sqrt(((i-i1)*(i-i1))+((i-i1)*(i-i1)))))/(SPAT_FILT_SIZE*SPAT_FILT_SIZE)); 
			 normgauss += (1/(sqrt(2*PI)*SPAT_FILT_SIZE)) * exp((-1/2)*((sqrt(((i-i1)*(i-i1))+((i-i1)*(i-i1))))*(sqrt(((i-i1)*(i-i1))+((i-i1)*(i-i1)))))/(SPAT_FILT_SIZE*SPAT_FILT_SIZE)); 
			 
			 /*}*/
		     }
		     
		     waveconvtmps[j][i]=waveconvtmps[j][i]/normgauss;		  
			
	   }
	}
	
        fp_grad=fopen(jac,"wb");

        /* output of the preconditioned gradient */
        for (i=1;i<=NXG;i++){
           for (j=1;j<=NYG;j++){
            
           fwrite(&waveconvtmps[j][i],sizeof(float),1,fp_grad);

           }
        }

        fclose(fp_grad);
	
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	fp_grad=fopen(jac,"rb");
	
	/* distribute spatial filtered gradient on computational nodes */
	for (i=1;i<=NXG;i++){
	   for (j=1;j<=NYG;j++){
	   
                        fread(&grad, sizeof(float), 1, fp_grad);
			   			
			if ((POS[1]==((i-1)/NX)) && 
		   	 (POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				waveconv[jj][ii]=grad;

			}
		}
	}
		
	fclose(fp_grad);

	free_matrix(waveconvtmp,1,NXG,1,NYG);
	free_matrix(waveconvtmps,1,NXG,1,NYG);

	
	} /* end of application condition for the spatial filter */


	
	
}



