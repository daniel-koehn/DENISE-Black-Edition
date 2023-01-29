/*
 * Apply spatial 2D Median filter to gradient/model
 *
 * Daniel Koehn
 * Kiel, 26/09/2016
 */

#include "fd.h"

void median_model(float ** waveconv, int filt_size)
{

	/* extern variables */

        extern float DH;
	extern int NX, NY, NXG, NYG, IDX, IDY;
	extern int MYID, MYID_SHOT, POS[3], FILT_SIZE;
	extern char JACOBIAN[STRING_SIZE];
	
	/* local variables */
	int i, j, ii, jj;
	int hfs;

	float **model_tmp, **kernel, **model_median, grad, smooth_meter;
	float sigma, s;
	
	char jac_tmp[STRING_SIZE];
	
	FILE *model, *FP1;
	
	char modfile[STRING_SIZE];

	
	/*if((MODEL_FILTER==1)||(GRAD_FILTER==1)){*/

		/* temporarily save gradient/model for 2D Median filtering */
	        sprintf(jac_tmp,"%s_median.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	        FP1=fopen(jac_tmp,"wb");
	                        
	        for (i=1;i<=NX;i=i+IDX){
	            for (j=1;j<=NY;j=j+IDY){
	                fwrite(&waveconv[j][i],sizeof(float),1,FP1);
	            }                       
	        }                           
                                    
	        fclose(FP1);                
                                    
	        MPI_Barrier(MPI_COMM_WORLD);
                                  
	        /* merge gradient/model file */ 
	        sprintf(jac_tmp,"%s_median.old",JACOBIAN);
	        if (MYID==0) mergemod(jac_tmp,3); 	

		MPI_Barrier(MPI_COMM_WORLD);

		if(MYID==0){			
			
			      	if (FILT_SIZE==0)	return;
			      	if (!(FILT_SIZE % 2)) {
			      	if (FILT_SIZE > 0)	FILT_SIZE += 1;
			      	else			FILT_SIZE -= 1;
			      	}				
	  
			      	hfs = abs(FILT_SIZE)/2;
				sigma = hfs/2;
				s = 2.0 * sigma * sigma;
			      	printf("\n hfs: %d \n",hfs);				
					      	
			      	model_tmp = matrix(-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
				model_median = matrix(1,NYG,1,NXG);
				kernel = matrix(1,abs(FILT_SIZE),1,abs(FILT_SIZE));
		      		      
				sprintf(jac_tmp,"%s_median.old",JACOBIAN);    
	
			      	model=fopen(jac_tmp,"rb");
			      	if (model==NULL) err(" Could not open gradient/model file !");
		
			      	/* load merged model */
			      	for (i=1;i<=NXG;i++){
					for (j=1;j<=NYG;j++){
						fread(&grad, sizeof(float), 1, model);
					      	model_tmp[j][i]=grad;
				      	}	
			      	}
		
			        fclose(model);
		      
			        /* apply 2D-Median filter*/
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
			
				/* 2D Median filter */
				for (j=1;j<=NYG;j++){
     				    for (i=1;i<=NXG;i++){

          			        /* write gridpoint data within the filter to kernel*/
	  				for (ii=-hfs;ii<=hfs;ii++){
	       			    		for (jj=-hfs;jj<=hfs;jj++){

							
							kernel[jj+hfs+1][ii+hfs+1] = model_tmp[j+jj][i+ii];

	               			    }
	          			}
					
					/* apply median filter */
					model_median[j][i] = median2d(kernel,abs(FILT_SIZE),abs(FILT_SIZE));																				

	      			    }
				}
			      
				/* output of median filterted gradients/models */      
				sprintf(jac_tmp,"%s_median.old",JACOBIAN);   
			      
				model=fopen(jac_tmp,"wb");
				for (i=1;i<=NXG;i++){
				  for (j=1;j<=NYG;j++){
			    
				  	fwrite(&model_median[j][i],sizeof(float),1,model);

				  }
				}
			
				fclose(model);
			
				free_matrix(model_tmp,-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
				free_matrix(model_median,1,NYG,1,NXG);
				free_matrix(kernel,1,abs(FILT_SIZE),1,abs(FILT_SIZE));
			
			} /* end of if(MYID==0)*/
		
		
		MPI_Barrier(MPI_COMM_WORLD);
		smooth_meter=FILT_SIZE*DH;
	
		if(MYID==0){printf("\n \t ---- 2D Median is applied to gradient/model (filter length of %d gridpoints which is equivalent to %4.2f meter) \n",FILT_SIZE,smooth_meter);}
	
		/* distribute smoothed jacobian on computational nodes */
		sprintf(jac_tmp,"%s_median.old",JACOBIAN);
		
		model=fopen(jac_tmp,"rb");
		if (model==NULL) err(" Could not open gradient/model file ! (distribute filtered gradient/model)");
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
	
		if(MYID==0){printf("\n \t ---- Filtered gradient/model is distributed on computational nodes ... ---- \n");}

        	/* clean up temporary files*/
        	MPI_Barrier(MPI_COMM_WORLD);
        	sprintf(jac_tmp,"%s_median.old.%i.%i",JACOBIAN,POS[1],POS[2]);
        	remove(jac_tmp);
	
	/*}*/

}
