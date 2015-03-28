/*------------------------------------------------------------------------
 *   Apply wavenumber domain filter to gradient
 *   
 *   Daniel Koehn
 *   Kiel, the 2nd of October 2012 
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include <fftw3.h>

void  wavenumber(float ** grad){

	/* declaration of extern variables */
        extern int NX, NY, NXG, NYG, IDX, IDY;
	extern int NPROCX, NPROCY, MYID, POS[3];
	extern int SPAT_FILT_SIZE;
	extern float WD_DAMP;
	extern char JACOBIAN[STRING_SIZE];
	
	/* declaration of local variables */
	int i,j, h, fa, fb, npadx, npady, itr;
	int tracl1, jj, ii, zeropad;
	float gradtmp;
	double npaddx, npaddy;
        float ** gradtmp1;
        
	char jac[STRING_SIZE];
	FILE *fp_grad, *FP3;
	                    
	/* size of the zeropadding layer */
        zeropad=SPAT_FILT_SIZE;
	                    
	/* temporarily save gradient for wavenumber filtering */
        sprintf(jac,"%s_wavenumber.old.%i%i",JACOBIAN,POS[1],POS[2]);
        FP3=fopen(jac,"wb");
                        
        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                fwrite(&grad[j][i],sizeof(float),1,FP3);
            }                       
        }                           
                                    
        fclose(FP3);                
                                    
        MPI_Barrier(MPI_COMM_WORLD);
                                  
        /* merge gradient file */ 
        sprintf(jac,"%s_wavenumber.old",JACOBIAN);
        if (MYID==0) mergemod(jac,3);  


if(MYID==0){    /* read the global model on node 0 and apply wavenumber damping */

        /*npadx = (int)(pow(2.0, ceil(log((double)(NXG))/log(2.0))+2.0) );
        npady = (int)(pow(2.0, ceil(log((double)(NYG))/log(2.0))+2.0) );*/
        
        npadx = NXG + 2.0*zeropad;
        npady = NYG + 2.0*zeropad;
        
        /* define temporary gradient matrix */
        gradtmp1 = matrix(1,npady,1,npadx);
        
        printf("npadx = %d \t npady = %d \n",npadx,npady);
        
	fftw_complex    *data, *fft_result, *ifft_result;
	fftw_plan       plan_forward, plan_backward;
                                        
	data        = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * npadx * npady );
	fft_result  = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * npadx * npady );
	ifft_result = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * npadx * npady );

        /* damping coefficient */
        /*damp=5e-5;
        damp=6e-5;*/

	printf("\n Spatial filter is applied to gradient (written by PE %d)\n",MYID); 
	
	fp_grad=fopen(jac,"rb");
	
	if (fp_grad==NULL) err(" Could not open gradient file ! ");
	
	/* load merged gradient */
	for (i=1;i<=npadx;i++){
	   for (j=1;j<=npady;j++){
	        
	        if((i>=1+zeropad)&&(i<=NXG+zeropad)&&(j>=1+zeropad)&&(j<=NYG+zeropad)){
	            fread(&gradtmp, sizeof(float), 1, fp_grad);
	            gradtmp1[j][i] = gradtmp;
		}
		else{
		    gradtmp1[j][i]=0.0;
		}
			
            }
	}
	
	fclose(fp_grad);
	
	/* Fill padding layer with non-zeros */
	/* top and bottom padding */
        for (i=1;i<=npadx;i++){
            for (j=1;j<=npady;j++){
            
                if((i>=1+zeropad)&&(i<=NXG+zeropad)&&(j<(1+zeropad))){
                  gradtmp1[j][i]=gradtmp1[1+zeropad][i];
                }
                
                if((i>=1+zeropad)&&(i<=NXG+zeropad)&&(j>(NYG+zeropad))){
                  gradtmp1[j][i]=gradtmp1[NYG+zeropad][i];
                }
            
            }
        }
        
        /* left and right padding */
        for (j=1;j<=npady;j++){
            for (i=1;i<=npadx;i++){
                                        
                if(i<(1+zeropad)){
                  gradtmp1[j][i]=gradtmp1[j][1+zeropad];
                }
                                                                                                          
                if(i>(NXG+zeropad)){                
                  gradtmp1[j][i]=gradtmp1[j][NXG+zeropad];         
                }
                                                                                                                                                                         
            }
        }
                                                                                                                                                                                             
	
        /* FFT of the  gradient */
        h=0;
        for (i=1;i<=npadx;i++){
            for (j=1;j<=npady;j++){
	                                                                                                               
                data[h][0] = gradtmp1[j][i];
                data[h][1] = 0.0;
	                                                                                                                                                                                                                     
                h++;
	                                                                                                                                                                                                                                                     
            }     
        }
	
	plan_forward  = fftw_plan_dft_2d(npady, npadx, data, fft_result, FFTW_FORWARD, FFTW_ESTIMATE );
        plan_backward = fftw_plan_dft_2d(npady, npadx, fft_result, ifft_result, FFTW_BACKWARD, FFTW_ESTIMATE );
	
        /* apply 2D-fft */
	fftw_execute( plan_forward );

	/* Apply Gaussian wavenumber damping */
        h=0;
        for (i=1;i<=npadx;i++){
            for (j=1;j<=npady;j++){ 
                
                if((i<=npadx/2+1)&&(j<=npady/2+1)){
                  fft_result[h][0] *=exp(-WD_DAMP*((i-1)*(i-1)+(j-1)*(j-1)));
                  fft_result[h][1] *=exp(-WD_DAMP*((i-1)*(i-1)+(j-1)*(j-1)));
                }
                
                if((i>npadx/2+1)&&(j<=npady/2+1)){ 
                  fft_result[h][0] *=exp(-WD_DAMP*((i-npadx)*(i-npadx)+(j-1)*(j-1)));
                  fft_result[h][1] *=exp(-WD_DAMP*((i-npadx)*(i-npadx)+(j-1)*(j-1)));
                }
                
                if((i<=npadx/2+1)&&(j>npady/2+1)&&(j<=npady)){
                  fft_result[h][0] *=exp(-WD_DAMP*((i-1)*(i-1)+(j-npady)*(j-npady)));          
                  fft_result[h][1] *=exp(-WD_DAMP*((i-1)*(i-1)+(j-npady)*(j-npady))); 
                }
                
                if((i>npadx/2+1)&&(j>npady/2+1)&&(j<=npady)){
                  fft_result[h][0] *=exp(-WD_DAMP*((i-npadx)*(i-npadx)+(j-npady)*(j-npady)));                            
                  fft_result[h][1] *=exp(-WD_DAMP*((i-npadx)*(i-npadx)+(j-npady)*(j-npady)));                  
                }
                                           
            h++;		   
	    }	   
        }
		
	/* apply 2D-ifft */
	fftw_execute( plan_backward );

	/* write damped gradient to temporary file */
	sprintf(jac,"%s_wavenumber.old",JACOBIAN);
	FP3=fopen(jac,"wb");
	h=0;
	for (i=1;i<=npadx;i++){
		for (j=1;j<=npady;j++){
			 
			 if((i>=1+zeropad)&&(i<=NXG+zeropad)&&(j>=1+zeropad)&&(j<=NYG+zeropad)){
			 gradtmp = ifft_result[h][0];
			 fwrite(&gradtmp,sizeof(float),1,FP3);}
			 
			 h++;
		}
	}
	fclose(FP3);
	        
	/* free memory */
	fftw_free( data );
	fftw_free( fft_result );
	fftw_free( ifft_result );

        fftw_destroy_plan( plan_forward );
        fftw_destroy_plan( plan_backward );

        free_matrix(gradtmp1,1,npady,1,npadx);
        
} /* end of if MYID==0*/

	 MPI_Barrier(MPI_COMM_WORLD);

         sprintf(jac,"%s_wavenumber.old",JACOBIAN);
	 FP3=fopen(jac,"rb");

	 /* distribute spatial filtered gradient on computational nodes */
	 for (i=1;i<=NXG;i++){
	    for (j=1;j<=NYG;j++){
			
			fread(&gradtmp, sizeof(float),1,FP3);

			if ((POS[1]==((i-1)/NX)) && 
		   	 (POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				grad[jj][ii]=gradtmp;

			}
			
		}
	}

        fclose(FP3);

        /* clean up temporary files*/
        MPI_Barrier(MPI_COMM_WORLD);
        sprintf(jac,"%s_wavenumber.old.%i%i",JACOBIAN,POS[1],POS[2]);
        remove(jac);
                                

}
