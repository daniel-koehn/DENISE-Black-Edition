/*
 * Apply spatial Gauss filter to gradient
 *
 * Daniel Koehn
 * Kiel, 26/09/2016
 */

#include "fd.h"

void gauss_filt(float ** waveconv)
{

	/* extern variables */

        extern float DH, WD_DAMP, FC_END, C_vs;
	extern int FREE_SURF, NX, NY, NXG, NYG, IDX, IDY;
	extern int NPROCX, NPROCY, MYID, POS[3];
	extern char JACOBIAN[STRING_SIZE];
	extern FILE *FP;
	extern int FILT_SIZE_GRAD;
	
	/* local variables */
	int i, j, ii, jj;
	int i1, j1, filtsize, hfs;

	float **model_tmp, **kernel, ** model_gauss, grad, normgauss, smooth_meter;
	float lam, sigma, s, sum, conv;
	
	char jac_tmp[STRING_SIZE];
	
	FILE *model, *FP1;
	
	char modfile[STRING_SIZE];

	/* temporarily save gradient for Gaussian filtering */
        sprintf(jac_tmp,"%s_gauss.old.%i%i",JACOBIAN,POS[1],POS[2]);
        FP1=fopen(jac_tmp,"wb");
                        
        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                fwrite(&waveconv[j][i],sizeof(float),1,FP1);
            }                       
        }                           
                                    
        fclose(FP1);                
                                    
        MPI_Barrier(MPI_COMM_WORLD);
                                  
        /* merge gradient file */ 
        sprintf(jac_tmp,"%s_gauss.old",JACOBIAN);
        if (MYID==0) mergemod(jac_tmp,3); 	

	if(MYID==0){
	
			lam = C_vs / FC_END;
			
			/* define filter size as fraction of reference velocity wavelength */
			FILT_SIZE_GRAD = round((WD_DAMP * lam)/DH);
			
		      	if (FILT_SIZE_GRAD==0)	return;
		      	if (!(FILT_SIZE_GRAD % 2)) {
		      	if (FILT_SIZE_GRAD > 0)	FILT_SIZE_GRAD += 1;
		      	else			FILT_SIZE_GRAD -= 1;
		      	}
	  
		      	hfs = abs(FILT_SIZE_GRAD)/2;
			sigma = hfs/2;
			s = 2.0 * sigma * sigma;
		      	printf("\n hfs: %d \n",hfs);
					      	
		      	model_tmp = matrix(-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
			model_gauss = matrix(-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
			kernel = matrix(1,abs(FILT_SIZE_GRAD),1,abs(FILT_SIZE_GRAD));
		      		      
			sprintf(jac_tmp,"%s_gauss.old",JACOBIAN);    
	
		      	model=fopen(jac_tmp,"rb");
		      	if (model==NULL) err(" Could not open gradient file !");
		
		      	/* load merged model */
		      	for (i=1;i<=NXG;i++){
				for (j=1;j<=NYG;j++){
					      fread(&grad, sizeof(float), 1, model);
				      	model_tmp[j][i]=grad;
			      	}	
		      	}
		
		        fclose(model);
		      
		      
		      
		        /* apply 2D-Gaussian filter on vp and vs model */
		        /* extrapolate array */
		      
		        /* left/right boundary */
		        for (j=1;j<=NYG;j++){
		      
		            for (i=-hfs+1;i<=0;i++){
			        model_tmp[j][i] = model_tmp[j][1];
			    }
			      
		            for (i=NXG+1;i<=NXG+hfs;i++){
			        model_tmp[j][i] = model_tmp[j][NXG];
			    }
			    
		        }
			      
		        /* top/bottom boundary incl. corners */
		        for (j=-hfs+1;j<=0;j++){
			
		             for (i=-hfs+1;i<=NXG+hfs;i++){
			         model_tmp[j][i] = model_tmp[1][i];
			     }
			     
		        }
		      
		        for (j=NYG+1;j<=NYG+hfs;j++){
			
		            for (i=-hfs+1;i<=NXG+hfs;i++){
			        model_tmp[j][i] = model_tmp[NYG][i];
			    }
			    
		        }
					
			/* create filter kernel */
			for (ii=-hfs;ii<=hfs;ii++){
			    for (jj=-hfs;jj<=hfs;jj++){
						      
			        kernel[jj+hfs+1][ii+hfs+1] = exp(-((ii*ii)/s) - ((jj*jj)/s));
				sum += kernel[jj+hfs+1][ii+hfs+1];
							      
			    }
			}
			
			/* normalize kernel */
			for (i=1;i<=FILT_SIZE_GRAD;i++){
     			    for (j=1;j<=FILT_SIZE_GRAD;j++){
         
         		        kernel[j][i] /= sum;

     			    }
			}
			
			/* apply Gaussian filter to gradient */
			for (j=1;j<=NYG;j++){
     			    for (i=1;i<=NXG;i++){

          		        conv = 0.0;
          		        /* loop over kernel*/
	  			for (ii=-hfs;ii<=hfs;ii++){
	       			    for (jj=-hfs;jj<=hfs;jj++){

	            		    conv += model_tmp[j+jj][i+ii] * kernel[jj+hfs+1][ii+hfs+1];				

               			    }
          			}

          			/* output of filtered gradient */
          			model_gauss[j][i] = conv;

      			    }
			}
			      
			/* output of smoothed gradients */      
			sprintf(jac_tmp,"%s_gauss.old",JACOBIAN);   
			      
			model=fopen(jac_tmp,"wb");
			for (i=1;i<=NXG;i++){
			  for (j=1;j<=NYG;j++){
			    
			  fwrite(&model_gauss[j][i],sizeof(float),1,model);

			  }
			}
			
			fclose(model);
			
			free_matrix(model_tmp,-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
			free_matrix(model_gauss,-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
			free_matrix(kernel,1,abs(FILT_SIZE_GRAD),1,abs(FILT_SIZE_GRAD));
			
		} /* end of if(MYID==0)*/
		
		
	MPI_Barrier(MPI_COMM_WORLD);
	smooth_meter=FILT_SIZE_GRAD*DH;
	
	if(MYID==0){printf("\n \t ---- Gradient is smoothed with Gaussian (filter length of %d gridpoints which is equivalent to %4.2f meter) \n",FILT_SIZE_GRAD,smooth_meter);}
	
	/* distribute smoothed jacobian on computational nodes */
	sprintf(jac_tmp,"%s_gauss.old",JACOBIAN);
		
	model=fopen(jac_tmp,"rb");
	if (model==NULL) err(" Could not open gradient file ! (distribute smoothed gradient)");
	for (i=1;i<=NXG;i++){
	   for (j=1;j<=NYG;j++){
	   
                        fread(&grad, sizeof(float), 1, model);
			   			
			if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				waveconv[jj][ii]=grad;

			}
		}
	}
		
	fclose(model);
	
	if(MYID==0){printf("\n \t ---- Smoothed gradient is distributed on computational nodes ... ---- \n");}

        /* clean up temporary files*/
        MPI_Barrier(MPI_COMM_WORLD);
        sprintf(jac_tmp,"%s_gauss.old.%i%i",JACOBIAN,POS[1],POS[2]);
        remove(jac_tmp);

}
