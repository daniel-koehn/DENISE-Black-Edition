/*
 * Apply spatial variable Gauss filter to gradient
 *
 * Daniel Koehn
 * Kiel, 18/09/2018
 */

#include "fd.h"

void gauss_filt_var(float ** waveconv, float ** vel_mod)
{

	/* extern variables */

        extern float DH, WD_DAMP, WD_DAMP1, FC_END, C_vs;
	extern int FREE_SURF, NX, NY, NXG, NYG, IDX, IDY;
	extern int NPROCX, NPROCY, MYID, POS[3];
	extern char JACOBIAN[STRING_SIZE];
	extern FILE *FP;
	extern int FILT_SIZE_GRAD;
	
	/* local variables */
	int i, j, ii, jj, FILT_SIZE_GRAD_X, FILT_SIZE_GRAD_Z;
	int i1, j1, filtsize, hfsx, hfsz, hfsmax;

	float **model_tmp, **grad_tmp, ** grad_gauss, grad, mod, normgauss;
	float lam, lam_max, sigmax, sx, sigmaz, sz, sum, conv, vel_max, kernel;
	
	char jac_tmp[STRING_SIZE];
	char modfile[STRING_SIZE];
	
	FILE *model, *FP1;
	
	/* temporarily save gradient for Gaussian filtering */
        sprintf(jac_tmp,"%s_gauss.old.%i.%i",JACOBIAN,POS[1],POS[2]);
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

	/* temporarily save velocity model for Gaussian filtering */
        sprintf(modfile,"%s_velmod.old.%i.%i",JACOBIAN,POS[1],POS[2]);
        FP1=fopen(modfile,"wb");
                        
        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                fwrite(&vel_mod[j][i],sizeof(float),1,FP1);
            }                       
        }                           
                                    
        fclose(FP1);                
                                    
        MPI_Barrier(MPI_COMM_WORLD);
                                  
        /* merge gradient file */ 
        sprintf(modfile,"%s_velmod.old",JACOBIAN);
        if (MYID==0) mergemod(modfile,3);

	if(MYID==0){

			/* load merged model */	
			model_tmp = matrix(1,NYG,1,NXG);      
			sprintf(modfile,"%s_velmod.old",JACOBIAN);    
	
		      	model=fopen(modfile,"rb");
		      	if (model==NULL) err(" Could not open model file !");
				      	
		      	for (i=1;i<=NXG;i++){
				for (j=1;j<=NYG;j++){
					fread(&mod, sizeof(float), 1, model);
				      	model_tmp[j][i] = mod;
			      	}	
		      	}
		
		        fclose(model);	

			/* calculate maximum velocity model value */
			vel_max = maximum_m(model_tmp, NXG, NYG);			
			
			/* define maximum filter size as fraction of maximum velocity wavelength */
			lam_max = vel_max / FC_END;
			FILT_SIZE_GRAD = round((WD_DAMP * lam_max)/DH);
			
		      	if (!(FILT_SIZE_GRAD % 2)) {
		      	    FILT_SIZE_GRAD += 1;
		      	}
			
			if (FILT_SIZE_GRAD < 6){FILT_SIZE_GRAD = 6;}
	  
		      	hfsmax = abs(FILT_SIZE_GRAD)/2;
			
			/* allocate temporary gradient matrices */		      	
		      	grad_tmp = matrix(-hfsmax+1,NYG+hfsmax,-hfsmax+1,NXG+hfsmax);
			grad_gauss = matrix(-hfsmax+1,NYG+hfsmax,-hfsmax+1,NXG+hfsmax);
		      	
			/* load merged gradient */	      
			sprintf(jac_tmp,"%s_gauss.old",JACOBIAN);    
	
		      	model=fopen(jac_tmp,"rb");
		      	if (model==NULL) err(" Could not open gradient file !");
				      	
		      	for (i=1;i<=NXG;i++){
				for (j=1;j<=NYG;j++){
					fread(&grad, sizeof(float), 1, model);
				      	grad_tmp[j][i]=grad;
			      	}	
		      	}
		
		        fclose(model);		      		      
		      
		        /* apply 2D-Gaussian filter on gradient */
		        /* extrapolate gradient array prior to filtering */
		      
		        /* left/right boundary */
		        for (j=1;j<=NYG;j++){
		      
		            for (i=-hfsmax+1;i<=0;i++){
			        grad_tmp[j][i] = grad_tmp[j][1];
			    }
			      
		            for (i=NXG+1;i<=NXG+hfsmax;i++){
			        grad_tmp[j][i] = grad_tmp[j][NXG];
			    }
			    
		        }
			      
		        /* top/bottom boundary incl. corners */
		        for (j=-hfsmax+1;j<=0;j++){
			
		             for (i=-hfsmax+1;i<=NXG+hfsmax;i++){
			         grad_tmp[j][i] = grad_tmp[1][i];
			     }
			     
		        }
		      
		        for (j=NYG+1;j<=NYG+hfsmax;j++){
			
		            for (i=-hfsmax+1;i<=NXG+hfsmax;i++){
			        grad_tmp[j][i] = grad_tmp[NYG][i];
			    }
			    
		        }							
			
			/* apply Gaussian filter to gradient */
			for (j=1;j<=NYG;j++){
     			    for (i=1;i<=NXG;i++){

				/* define maximum filter size as fraction of local velocity wavelength */
				lam = model_tmp[j][i] / FC_END;
				
				FILT_SIZE_GRAD_X = round((WD_DAMP * lam)/DH);
				FILT_SIZE_GRAD_Z = round((WD_DAMP1 * lam)/DH);
			
		      		if (FILT_SIZE_GRAD_Z>0){
				
		      		    if (!(FILT_SIZE_GRAD_X % 2)) {
		      		        FILT_SIZE_GRAD_X += 1;
		      		    }

		      		    if (!(FILT_SIZE_GRAD_Z % 2)) {
		      		        FILT_SIZE_GRAD_Z += 1;
		      		    }
				    
				    if (FILT_SIZE_GRAD_X < 6){FILT_SIZE_GRAD_X = 6;}
	  		            if (FILT_SIZE_GRAD_Z < 6){FILT_SIZE_GRAD_Z = 6;}
	  
		      		    hfsx = abs(FILT_SIZE_GRAD_X)/2;
				    sigmax = hfsx/2;
				    sx = 2.0 * sigmax * sigmax;

		      		    hfsz = abs(FILT_SIZE_GRAD_Z)/2;
				    sigmaz = hfsz/2;
				    sz = 2.0 * sigmaz * sigmaz;

          		            conv = 0.0;
				    sum = 0.0;
          		            /* loop over kernel*/
	  			    for (ii=-hfsx;ii<=hfsx;ii++){
	       			        for (jj=-hfsz;jj<=hfsz;jj++){

					    kernel = exp(-((ii*ii)/sx) - ((jj*jj)/sz)); 
	            		    	    conv += grad_tmp[j+jj][i+ii] * kernel;
					    sum += kernel;				

               			        }
          			    }

          			    /* output of filtered gradient */
          			    grad_gauss[j][i] = conv/sum;				    
				    
			       }else{
			      
			            grad_gauss[j][i] = grad_tmp[j][i];
			      
			       }			        	    

      			    }
			}
			      
			/* output of smoothed gradients */      
			sprintf(jac_tmp,"%s_gauss.old",JACOBIAN);   
			      
			model=fopen(jac_tmp,"wb");
			for (i=1;i<=NXG;i++){
			  for (j=1;j<=NYG;j++){
			    
			      fwrite(&grad_gauss[j][i],sizeof(float),1,model);

			  }
			}
			
			fclose(model);
			
			free_matrix(model_tmp,1,NYG,1,NXG);
			free_matrix(grad_tmp,-hfsmax+1,NYG+hfsmax,-hfsmax+1,NXG+hfsmax);
			free_matrix(grad_gauss,-hfsmax+1,NYG+hfsmax,-hfsmax+1,NXG+hfsmax);
			
		} /* end of if(MYID==0)*/
		
		
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(MYID==0){printf("\n \t ---- Gradient is smoothed with spatial variable Gaussian filter\n");}
	
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
        sprintf(jac_tmp,"%s_gauss.old.%i.%i",JACOBIAN,POS[1],POS[2]);
        remove(jac_tmp);

}
