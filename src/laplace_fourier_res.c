/*-----------------------------------------------------------------------------
 *   Calculate data residuals and objective function in Laplace-Fourier domain 
 *
 *   Daniel Koehn
 *   Kiel, the 14th of May 2015
 *
 *  ---------------------------------------------------------------------------*/

#include "fd.h"
#include <complex.h>
#include <fftw3.h>

void laplace_fourier_res(float **sectionvy_obs, float **sectionvy, float **sectiondiff, int ntr, int ntr_glob, int ns, int ishot, int nshots, int iter, int **recpos, int **recpos_loc, float **srcpos){ 

	/* declaration of global variables */
	extern float DT, DH, OFFSETC;
	extern int SEIS_FORMAT, MYID, NT, TIMEWIN;
	extern char  SEIS_FILE_VY[STRING_SIZE], PARA[STRING_SIZE], DATA_DIR[STRING_SIZE];
	extern int TRKILL, OFFSET_MUTE, GRAD_FORM;
	extern char TRKILL_FILE[STRING_SIZE];
	
	/* declaration of variables for trace killing */
	int ** kill_tmp, *kill_vector, h, j, i, Npad;
	char trace_kill_file[STRING_SIZE];	
	FILE *ftracekill;

        /* declaration of variables for offset-muting */
        float offset, xr, yr, xs, ys;

        /* complex variables for source wavelet estimation */
        fftw_complex *D_ss, *D_ss_fd;

        /* parameters for STA/LTA first arrival picking */
        float *picked_times=NULL;
	
	/* variables for data integration */
	int invtime;
	float **integrated_section=NULL, **integrated_sectiondata=NULL;
	float EPS_NORM;
	
        integrated_section = matrix(1,ntr,1,ns);
	integrated_sectiondata = matrix(1,ntr,1,ns);
	
	Npad = 2.0 * ns;
        Npad = (int)(pow(2.0, ceil(log((double)(Npad))/log(2.0))+2.0) );
    
        EPS_NORM=1e0;
    
        /* Allocate memory */
        D_ss  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
        D_ss_fd  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
        
        
	/* reverse time direction integrate data if GRAD_FORM=1 */
        for(i=1;i<=ntr;i++){
    
           invtime = ns;
	   for(j=1;j<=ns;j++){
	   
	   
	      if(GRAD_FORM==1){
	         integrated_section[i][invtime] += DT*sectionvy[i][j];
	         integrated_sectiondata[i][invtime] += DT*sectionvy_obs[i][j];
	      }
	      
	      if((GRAD_FORM==2)||(GRAD_FORM==3)){
	         integrated_section[i][invtime] = sectionvy[i][j];
	         integrated_sectiondata[i][invtime] = sectionvy_obs[i][j];
	      }	      
	      
	      invtime--;
	   }	
	   
	}
		
        /* pick first arrivals in the synthetic data with STA/LTA-picker and apply time window to field data before stf inversion */ 	                                           
        /*if(TIMEWIN==3){
	  picked_times = vector(1,ntr);
	  stalta(sectionvy, ntr, ns, picked_times, ishot);
          time_window_stf(sectionvy_obs, iter, ntr, ns, ishot);
        }*/
	
	TIMEWIN=2;
	
	/* apply time damping to field and model data */
	picked_times = vector(1,ntr);
	time_window(integrated_sectiondata,picked_times,iter,ntr_glob,recpos_loc,ntr,ns,ishot);
        time_window(integrated_section,picked_times,iter,ntr_glob,recpos_loc,ntr,ns,ishot);
	                  
			       
	/* TRKILL==1 - trace killing is applied */
	if(TRKILL){
	  kill_tmp = imatrix(1,nshots,1,ntr);
	  kill_vector = ivector(1,ntr);

	  ftracekill=fopen(TRKILL_FILE,"r");

	  if (ftracekill==NULL) err(" Trace kill file could not be opened!");

		for(i=1;i<=nshots;i++){
			for(j=1;j<=ntr;j++){
				fscanf(ftracekill,"%d",&kill_tmp[i][j]);
			}
		}

		fclose(ftracekill);

		for(i=1;i<=ntr;i++){
	   	   kill_vector[i] = kill_tmp[ishot][i];
		}

	  for(i=1;i<=ntr;i++){

	     if(kill_vector[i]==1){
	       
	        for(j=1;j<=ns;j++){
		   integrated_section[i][j]=0.0;
		   integrated_sectiondata[i][j]=0.0;
	        }
	     }	
    	     
	  }
	
	}
		
	/* trace killing ends here */

        /* apply offset mute */
        if(OFFSET_MUTE){

         /*printf("OFFSETC = %f \n",OFFSETC);
         printf("OFFSET_MUTE = %d \n",OFFSET_MUTE);      */

         for (i=1;i<=ntr;i++){
             
             /* calculate source and receiver positions */
      	     xr = recpos[1][recpos_loc[3][i]]*DH;
      	     xs = srcpos[1][ishot];
      	     yr = recpos[2][recpos_loc[3][i]]*DH;
      	     ys = srcpos[2][ishot];

             /* calculate absolute offset */
             offset = sqrt(((xs-xr)*(xs-xr))+((ys-yr)*(ys-yr)));

             /* mute far-offset data*/
             if((OFFSET_MUTE==1)&&(offset>=OFFSETC)){
                
                for(j=1;j<=ns;j++){
		   integrated_section[i][j]=0.0;
		   integrated_sectiondata[i][j]=0.0;
	        }
                    
             }

             /* mute near-offset data*/
             if((OFFSET_MUTE==2)&&(offset<=OFFSETC)){

	        for(j=1;j<=ns;j++){
		   integrated_section[i][j]=0.0;
		   integrated_sectiondata[i][j]=0.0;
	        }

             }
         } 

        } /* end of OFFSET_MUTE */	                                 


        /* FFT of each data and model trace and calculation of nominator and denominator of Wiener deconvolution */
        for(i=1;i<=ntr;i++){

           /* allocate memory for complex variables */           
           fftw_complex *in_data, *out_data, *in_model, *out_model, *res_td, *res;
           fftw_plan p_data,p_model;
         
           in_data  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
           out_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);

           in_model  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
           out_model = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
	   
	   res_td  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
	   res  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);

           /* define real and imaginary parts of data vectors and apply zero-padding */
           for(j=0;j<Npad;j++){

              if(j<ns){
                 
	      	 in_model[j] = (integrated_section[i][j+1]) + 0.0*I;
                 in_data[j]  = (integrated_sectiondata[i][j+1]) + 0.0*I;	
              }
              else{
	      	in_model[j] = 0.0 + 0.0*I;
                 in_data[j] = 0.0 + 0.0*I;
              }
              
	   }
           
           /* apply FFTW */           
           p_data  = fftw_plan_dft_1d(Npad, in_data, out_data, 1, FFTW_ESTIMATE);
           p_model = fftw_plan_dft_1d(Npad, in_model, out_model, 1, FFTW_ESTIMATE);         

           fftw_execute(p_data);
           fftw_execute(p_model);

           /* calculate data residuals in the Laplace-Fourier domain (Jun etal., 2014)*/
           for(j=0;j<Npad;j++){

              /* real parts of the nominator and denominator */
              res[j] = (1.0/out_model[j]+EPS_NORM)*conj(log(out_model[j]+EPS_NORM)-log(out_data[j]+EPS_NORM));
	      
           }
	   
	   /* inverse FFTW of data residuals */
           fftw_plan p_stf;
           p_stf  = fftw_plan_dft_1d(Npad, res, res_td, -1, FFTW_ESTIMATE);
           fftw_execute(p_stf);

           /* extract real part and write data to sectiondiff */
	   /*h=1;
           for(j=Npad;j>Npad-ns;j--){
       	      sectiondiff[i][h]=creal(res_td[h])/Npad;
	      h=h+1;
           }*/  
	   
	 
	   for(j=1;j<=ns;j++){
	      sectiondiff[i][j] = integrated_section[i][j] - integrated_sectiondata[i][j];    	
	   }         

           fftw_destroy_plan(p_data);
           fftw_free(in_data); 
           fftw_free(out_data);          

           fftw_destroy_plan(p_model);
	   fftw_destroy_plan(p_stf);
           fftw_free(in_model); 
           fftw_free(out_model);
	   fftw_free(res);
	   fftw_free(res_td);
	   
        }
        						
	/* free memory for trace killing and FFTW */
	if(TRKILL){
	   free_imatrix(kill_tmp,1,nshots,1,ntr);
	   free_ivector(kill_vector,1,ntr);
	}

         
 	fftw_free(D_ss);
	fftw_free(D_ss_fd);
	
	
	free_vector(picked_times,1,ntr);
	free_matrix(integrated_section,1,ntr,1,ns);
        free_matrix(integrated_sectiondata,1,ntr,1,ns);
 
}

