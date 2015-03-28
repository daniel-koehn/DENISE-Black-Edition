/*
 * taper gradient with a gaussian frame to damp inversion artefacts near the sources and receivers
 sws == 1 vertical taper (for tomography geometry)
 sws == 2 horizontal taper (for reflection geometry)
 sws == 3 local radial taper at the source and receiver positions
 sws == 4  
 */

#include "fd.h"

void spat_agc(float ** waveconv, int iter, int sws)
{

	/* extern variables */

        extern float DH;
	extern int FREE_SURF, NX, NY, NXG, NYG;
	extern int NPROCX, NPROCY, MYID, POS[3];
	extern char JACOBIAN[STRING_SIZE];
	extern FILE *FP;
	
	/* local variables */
	int i, j, h, ifw, ii, jj, n, xb, yb, xe, ye, taperlength,taperlength2, SPAT_FILT_SIZE, SPAT_FILT_1;
	int ijc, iy, ix, iii, jjj, xx, yy, srctaper_gridpt, i1, j1, filtsize;

	float **waveconvtmp, **waveconvtmps, grad, normgauss;
	int winx, winy;
	float min2, max2, min1, max1, norm;
	const float eps=1e-20; /* Levenberg parameter */
	char jac[STRING_SIZE];
	FILE *fp_grad;
	
	SPAT_FILT_SIZE=5;
	winx = SPAT_FILT_SIZE;
	winy = SPAT_FILT_SIZE;
	
        waveconvtmp = matrix(1,NYG,1,NXG);
	waveconvtmps = matrix(1,NYG,1,NXG);
	
	if(sws==1){
	sprintf(jac,"%s_g.old",JACOBIAN);}

	if(sws==2){
	sprintf(jac,"%s_g_u.old",JACOBIAN);}
	
	if(sws==3){
	sprintf(jac,"%s_g_rho.old",JACOBIAN);}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(MYID==0){
	fp_grad=fopen(jac,"rb");
	if (fp_grad==NULL) err(" Could not open gradient file ! ");
	
	/* load merged gradient */
	for (i=1;i<=NXG;i++){
	   for (j=1;j<=NYG;j++){
                  fread(&grad, sizeof(float), 1, fp_grad);
			waveconvtmp[j][i]=grad;
			waveconvtmps[j][i]=0.0;
		}
	}
	
	fclose(fp_grad);
	
	min2 = waveconvtmp[1][1];
	max2 = waveconvtmp[1][1];
	
	/* estimate min and max of the global gradient */
	 for (i=1;i<=NXG;i++){
	            for (j=1;j<=NYG;j++){
	            
	             if(waveconvtmp[j][i]<min2){min2=waveconvtmp[j][i];}
	             if(waveconvtmp[j][i]>max2){max2=waveconvtmp[j][i];}
	            
	            }
	 }
	
	/* Normalize gradient */
	norm = sqrt(max2*max2);
	
	printf("Norm = %e \n",norm);
	
	if(sqrt(min2*min2)>sqrt(max2*max2)){norm=sqrt(min2*min2);}
	
	for (i=1;i<=NXG;i++){
	    for (j=1;j<=NYG;j++){	    	    
	       waveconvtmp[j][i] = waveconvtmp[j][i]/norm;	    
	    }
	}                
	
	/* estimate Min and Max in AGC window */
	for (i=1;i<=NXG;i=i+winx){
	   for (j=1;j<=NYG;j=j+winy){

                     max1 = min2;   
		     min1 = max2;
		     
		     for (ii=i-winx;ii<=i+winx;ii++){
	                 for (jj=j-winy;jj<=j+winy;jj++){
			 
	                     if((jj>0)&&(jj<NYG+1)&&(ii>0)&&(ii<NXG+1)){
	                      
	                        if(waveconvtmp[jj][ii]>max1){max1=waveconvtmp[jj][ii];}
	                        if(waveconvtmp[jj][ii]<min1){min1=waveconvtmp[jj][ii];}
	                        
			     }
			  }
		     }

                 if(sqrt(max1*max1)<sqrt(min1*min1)){max1=min1;}		     

 		     for (ii=i-winx;ii<=i+winx;ii++){
	                 for (jj=j-winy;jj<=j+winy;jj++){
			 
	                        if((jj>0)&&(jj<NYG+1)&&(ii>0)&&(ii<NXG+1)){
	                        waveconvtmps[jj][ii] = waveconvtmp[jj][ii]/sqrt((max1+eps)*(max1+eps));
	                        }
			 }
		     }
                 
				
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


	
	
}



