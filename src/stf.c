/*------------------------------------------------------------------------
 *   Inversion of source time function 
 *
 *   Daniel Koehn
 *   Kiel, the 12th of September 2013
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include <complex.h>
#include <fftw3.h>

void stf(float **sectionvy_obs, float **sectionvy, int ntr_glob, int ishot, int ns, int iter, int nshots, float **signals, int **recpos, float **srcpos){ 

	/* declaration of global variables */
	extern float DT, DH, OFFSETC, EPS_STF, OFFSETC_STF;
	extern int SEIS_FORMAT, MYID, NT, NORMALIZE, TIMEWIN, INV_STF;
	extern char  SEIS_FILE_VY[STRING_SIZE], PARA[STRING_SIZE], DATA_DIR[STRING_SIZE];
	extern int TRKILL, OFFSET_MUTE;
	extern char TRKILL_FILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE];
	
	/* declaration of variables for trace killing */
        float *STF_vector;
	int ** kill_tmp, *kill_vector, h, j, i, Npad;
	char trace_kill_file[STRING_SIZE];	
        double npadd;
	FILE *ftracekill, *STF;

        /* declaration of variables for offset-muting */
        float offset, xr, yr, xs, ys;

        /* complex variables for source wavelet estimation */
        fftw_complex *sumn, *sumd, *D_s, *D_ss, *D_ss_fd, *D_s_td;
        float Ebar;
        char signal_wave[STRING_SIZE];

        /* parameters for STA/LTA first arrival picking */
        float *picked_times=NULL;

        sprintf(signal_wave,"%s_shot_%i.dat",SIGNAL_FILE,ishot);
	
	printf("\n================================================================================================\n\n");
	printf("\n ***** Inversion of Source Time Function - shot: %d - it: %d ***** \n\n",ishot,iter);

        Npad = (int)(pow(2.0, ceil(log((double)(ns))/log(2.0))+2.0) );
        /*Npad = ns;*/
    
        /* Allocate memory */
        sumn  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
        sumd  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
        D_s  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
        D_ss  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
        D_ss_fd  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
        D_s_td  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);

        STF_vector = vector(1,ns);

        /*printf("Npad=%d \n",Npad);   */
        
        /* pick first arrivals in the synthetic data with STA/LTA-picker and apply time window to field data before stf inversion */ 	                                           
        if(INV_STF==2){
	  picked_times = vector(1,ntr_glob);
	  stalta(sectionvy, ntr_glob, ns, picked_times, ishot);
          time_window_stf(sectionvy_obs, iter, ntr_glob, ns, ishot);
        }
                               
	/* TRKILL==1 - trace killing is applied */
	if(TRKILL){
	  kill_tmp = imatrix(1,nshots,1,ntr_glob);
	  kill_vector = ivector(1,ntr_glob);

	  ftracekill=fopen(TRKILL_FILE,"r");

	  if (ftracekill==NULL) err(" Trace kill file could not be opened!");

		for(i=1;i<=nshots;i++){
			for(j=1;j<=ntr_glob;j++){
				fscanf(ftracekill,"%d",&kill_tmp[i][j]);
			}
		}

		fclose(ftracekill);

		for(i=1;i<=ntr_glob;i++){
	   	   kill_vector[i] = kill_tmp[ishot][i];
		}

	}
	
	if(TRKILL){
	  for(i=1;i<=ntr_glob;i++){

	     if(kill_vector[i]==1){
	       
	        for(j=1;j<=ns;j++){
		   sectionvy[i][j]=0.0;
		   sectionvy_obs[i][j]=0.0;
	        }
	     }	
    	     
	  }
	
	}	
	/* trace killing ends here */

        /* apply offset mute */
        if(OFFSET_MUTE){

         /*printf("OFFSETC = %f \n",OFFSETC);
         printf("OFFSET_MUTE = %d \n",OFFSET_MUTE);      */

         for (i=1;i<=ntr_glob;i++){
             
             /* calculate source and receiver positions */
      	     xr = recpos[1][i]*DH;
      	     xs = srcpos[1][ishot];
      	     yr = recpos[2][i]*DH;
      	     ys = srcpos[2][ishot];

             /* calculate absolute offset */
             offset = sqrt(((xs-xr)*(xs-xr))+((ys-yr)*(ys-yr)));

             printf("offset = %f \n",offset);

             /* mute far-offset data*/
             if((OFFSET_MUTE==1)&&(offset>=OFFSETC)){
                
                for(j=1;j<=ns;j++){
		   sectionvy[i][j]=0.0;
		   sectionvy_obs[i][j]=0.0;
	        }
                    
             }

             /* mute near-offset data*/
             if((OFFSET_MUTE==2)&&(offset<=OFFSETC)){

	        for(j=1;j<=ns;j++){
		   sectionvy[i][j]=0.0;
		   sectionvy_obs[i][j]=0.0;
	        }

             }
         } 

        } /* end of OFFSET_MUTE */	

        /* apply offset mute for STF inversion */
        if(OFFSETC_STF>0.0){

         for (i=1;i<=ntr_glob;i++){
             
             /* calculate source and receiver positions */
      	     xr = recpos[1][i]*DH;
      	     xs = srcpos[1][ishot];
      	     yr = recpos[2][i]*DH;
      	     ys = srcpos[2][ishot];

             /* calculate absolute offset */
             offset = sqrt(((xs-xr)*(xs-xr))+((ys-yr)*(ys-yr)));

             /* mute far-offset data*/
             if(fabs(offset)>=OFFSETC_STF){

	        for(j=1;j<=ns;j++){
		   sectionvy[i][j]=0.0;
		   sectionvy_obs[i][j]=0.0;
	        }

             }
         } 

        } /* end of OFFSETC_STF */

        /* Trace normalization to maximum amplitude of each trace */
        if(NORMALIZE==1){
          normalize_data(sectionvy,ntr_glob,ns);
          normalize_data(sectionvy_obs,ntr_glob,ns);
        }
        
        /* Trace normalization of field data with respect to maximum amplitude of model data */
        if(NORMALIZE==2){
          normalize_data_rel(sectionvy,sectionvy_obs,ntr_glob,ns);
        }
                                          
        /* initialize nominator and denominator for Wiener deconvolution */	
        for(j=0;j<Npad;j++){

           sumn[j]=0.0 + 0.0*I;
           sumd[j]=0.0 + 0.0*I;

        }

        /* FFT spike wavelet */
        for(j=0;j<Npad;j++){

           if(j<ns){
	     D_ss[j] = (double) (signals[1][j+1]) + 0.0*I;
           }
           else{D_ss[j] = 0.0 + 0.0*I;}          

        }

        fftw_plan p_s;

        p_s = fftw_plan_dft_1d(Npad, D_ss, D_ss_fd, 1, FFTW_ESTIMATE);         
        fftw_execute(p_s);

        /* FFT of each data and model trace and calculation of nominator and denominator of Wiener deconvolution */
        Ebar = 0.0;
        for(i=1;i<=ntr_glob;i++){

           /* allocate memory for complex variables */           
           fftw_complex *in_data, *out_data, *in_model, *out_model;
           fftw_plan p_data,p_model;
         
           in_data  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
           out_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);

           in_model  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
           out_model = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);

           /* define real and imaginary parts of data vectors and apply zero-padding */
           for(j=0;j<Npad;j++){

              if(j<ns){
                 
	      	 in_model[j] = (double) (sectionvy[i][j+1]) + 0.0*I;
                 in_data[j]  = (double) (sectionvy_obs[i][j+1]) + 0.0*I;	
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

           /* estimate nominator and denominator of the Wiener deconvolution */
           for(j=0;j<Npad;j++){

              /* real parts of the nominator and denominator */
              sumn[j] += ntr_glob*out_data[j]*conj(out_model[j]);
              sumd[j] += ntr_glob*out_model[j]*conj(out_model[j]);
              Ebar += sumd[j];
           }           

           fftw_destroy_plan(p_data);
           fftw_free(in_data); 
           fftw_free(out_data);          

           fftw_destroy_plan(p_model);
           fftw_free(in_model); 
           fftw_free(out_model);

        }

        Ebar = Ebar/(Npad*ntr_glob);

        /* construct source wavelet in frequency domain by Wiener deconvolution */        		
        for(j=0;j<Npad;j++){
           
          D_s[j] = D_ss_fd[j]*(sumn[j]/(sumd[j]+(EPS_STF*ntr_glob*Ebar)));
          /*D_s[j] = (sumn[j]/(sumd[j]+(EPS_STF*ntr_glob*Ebar)));*/          
            
        }        

        /* inverse FFTW of the estimated STF */
        fftw_plan p_stf;
        p_stf  = fftw_plan_dft_1d(Npad, D_s, D_s_td, -1, FFTW_ESTIMATE);
        fftw_execute(p_stf);

        /* extract real part and flip STF_vector */
        for(j=0;j<ns;j++){
           STF_vector[j+1]=creal(D_s_td[j])/Npad;
        }

        /* FFT of STF */
        /*for(j=0;j<Npad;j++){

           if(j<ns){
	     D_s_td[j] = (double) (STF_vector[j+1]);
           }
           else{D_ss[j] = 0.0;}          

        }

        fftw_plan STF_s;

        STF_s = fftw_plan_dft_1d(Npad, D_s_td, D_s, 1, FFTW_ESTIMATE);         
        fftw_execute(STF_s);*/

        /* convolve STF with spike */
        /*for(j=0;j<Npad;j++){

           D_s[j] = D_ss_fd[j]*D_s[j];
                      
        }*/

        /* inverse FFTW of the estimated STF */
        /*p_stf  = fftw_plan_dft_1d(Npad, D_s, D_s_td, -1, FFTW_ESTIMATE);
        fftw_execute(p_stf);*/

        /* extract real part of STF*Spike vector */
        /*for(j=0;j<ns;j++){
           STF_vector[j+1]=creal(D_s_td[j+1]);
        }*/
	
        /* normalization of source wavelet to maximum amplitude for better numerical performance */
        /*normalize_STF(STF_vector,ns);*/

        /* output of the STF */
        STF=fopen(signal_wave,"w");                                	
        for(j=1;j<=ns;j++){
           fprintf(STF,"%e\n",STF_vector[j]);
        }	
        fclose(STF);
						
	printf("\n\n================================================================================================\n");
	
	
	/* free memory for trace killing and FFTW */
	if(TRKILL){
	free_imatrix(kill_tmp,1,nshots,1,ntr_glob);
	free_ivector(kill_vector,1,ntr_glob);
	}

        fftw_free(sumn); 
        fftw_free(sumd); 
	fftw_free(D_s);
 	fftw_free(D_ss);
	fftw_free(D_ss_fd);
	fftw_free(D_s_td);
        fftw_destroy_plan(p_stf);
        fftw_destroy_plan(p_s);
	
        /*fftw_destroy_plan(STF_s);*/

        free_vector(STF_vector,1,ns);
	
	if(INV_STF==2){
	  free_vector(picked_times,1,ntr_glob);
	}  

}

