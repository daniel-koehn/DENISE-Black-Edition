/*------------------------------------------------------------------------
 *   Calculate Data Residuals                                  
 *   last update 29/03/08, D.Koehn
 *  ----------------------------------------------------------------------*/
#include "fd.h"

double calc_res(float **sectiondata, float **section, float **sectiondiff, float **sectiondiffold, int ntr, int ns, int LNORM, float L2, int itest, int sws, int swstestshot, int ntr_glob, int **recpos, int **recpos_loc, float **srcpos, int nsrc_glob, int ishot, int iter){

/* declaration of variables */
extern float DT, DH, OFFSETC, FC, FC_START;
extern int REC1, REC2, MYID, ORDER, INVMAT;
extern int TRKILL, TIME_FILT, GRAD_FORM,ENV;
extern char TRKILL_FILE[STRING_SIZE];
extern int NORMALIZE, TIMEWIN, RTM, OFFSET_MUTE;
float RMS, RMS_obs, signL1, intseis;
int Lcount,i,j,invtime,k,h, umax=0;
float l2;
float abs_data, abs_synthetics, data_mult_synthetics, intseis_data, intseis_synthetics;
float intseis_section, intseis_sectiondata, offset, xr, yr, xs, ys;
float *picked_times=NULL, eps;
float **integrated_section=NULL, **integrated_sectiondata=NULL;
float **integrated_sectiondata_envelope=NULL, **integrated_section_envelope=NULL, **integrated_section_hilbert=NULL;
float **dummy_1=NULL, **dummy_2=NULL; 
float EPS_LNORM, EPS_LNORM6, tmp, tmp1;

/* EPS_LNORM=1e1; */
EPS_LNORM=0.0;
EPS_LNORM6=1e-1;
eps=1e-15;

if(TIMEWIN) picked_times = vector(1,ntr);

/* sectiondiff will be set to zero*/
umax=ntr*ns;
zero(&sectiondiff[1][1],umax);


/* declaration of variables for trace killing */
int ** kill_tmp, *kill_vector;
char trace_kill_file[STRING_SIZE];	
FILE *ftracekill;
if(TRKILL){
  
	/* sectiondiff will be set to zero*/
	umax=ntr*ns;
	zero(&sectiondiff[1][1],umax);
	
	kill_tmp = imatrix(1,nsrc_glob,1,ntr_glob);
	kill_vector = ivector(1,ntr);

	ftracekill=fopen(TRKILL_FILE,"r");

	if (ftracekill==NULL) err(" Trace kill file could not be opened!");

	for(i=1;i<=nsrc_glob;i++){
		for(j=1;j<=ntr_glob;j++){
			fscanf(ftracekill,"%d",&kill_tmp[i][j]);
		}
	}

	fclose(ftracekill);

	h=1;
	for(i=1;i<=ntr;i++){
	   kill_vector[h] = kill_tmp[ishot][recpos_loc[3][i]];
	   h++;
	}
} /* end if(TRKILL)*/


RMS=0.0;
Lcount=1;  

if((LNORM==5)&&(GRAD_FORM==1)){

	integrated_section = matrix(1,ntr,1,ns);
	integrated_sectiondata = matrix(1,ntr,1,ns);

        /* Integration of measured and synthetic data  */
        for(i=1;i<=ntr;i++){
                intseis_section = 0.0;
      	        for(j=1;j<=ns;j++){
		        intseis_section += section[i][j];
		        integrated_section[i][j]=intseis_section;
		}
	}

        for(i=1;i<=ntr;i++){
               intseis_sectiondata = 0.0;
      	       for(j=1;j<=ns;j++){
		        intseis_sectiondata += sectiondata[i][j];
		        integrated_sectiondata[i][j]=intseis_sectiondata;
		}
	}

        /*if((TIMEWIN==1)||(TIMEWIN==2)){
            time_window(integrated_sectiondata, picked_times, iter, ntr_glob,recpos_loc, ntr, ns, ishot);
            time_window(integrated_section, picked_times, iter, ntr_glob,recpos_loc, ntr, ns, ishot);
        }*/

} /* end of if LNORM==5 */

/* Envelope objective functions after Chi, Dong and Liu (2014) */
if(LNORM==6){

	integrated_section = matrix(1,ntr,1,ns);
	integrated_sectiondata = matrix(1,ntr,1,ns);

        integrated_section_envelope = matrix(1,ntr,1,ns); 
        integrated_sectiondata_envelope = matrix(1,ntr,1,ns); 
        integrated_section_hilbert = matrix(1,ntr,1,ns); 
        dummy_1 = matrix(1,ntr,1,ns);    
        dummy_2 = matrix(1,ntr,1,ns); 


        /* Integration of measured and synthetic data  */
        for(i=1;i<=ntr;i++){
                intseis_section = 0.0;
      	        for(j=1;j<=ns;j++){
		
		        if(GRAD_FORM==1){
		          intseis_section += section[i][j];
		          integrated_section[i][j]=intseis_section;
			}
			
		        if((GRAD_FORM==2)||(GRAD_FORM==3)){
		          integrated_section[i][j]=section[i][j];
			}
			
		}
	}

        for(i=1;i<=ntr;i++){
               intseis_sectiondata = 0.0;
      	       for(j=1;j<=ns;j++){
	       
	       		if(GRAD_FORM==1){
		          intseis_sectiondata += sectiondata[i][j];
		          integrated_sectiondata[i][j]=intseis_sectiondata;
			}
			
		        if((GRAD_FORM==2)||(GRAD_FORM==3)||(GRAD_FORM==4)){
		          integrated_sectiondata[i][j]=sectiondata[i][j];
			}
	       
		}
	}
	
	calc_envelope(integrated_sectiondata,integrated_sectiondata_envelope,ns,ntr); 
        calc_envelope(integrated_section,integrated_section_envelope,ns,ntr);  
        calc_hilbert(integrated_section,integrated_section_hilbert,ns,ntr); 

	/* low-pass filter envelope data */
	timedomain_filt(integrated_sectiondata_envelope,FC,ORDER,ntr,ns,1);
	timedomain_filt(integrated_section_envelope,FC,ORDER,ntr,ns,1);
	timedomain_filt(integrated_section_hilbert,FC,ORDER,ntr,ns,1);


        /* L2 envelope objective function*/
        if(ENV==1){

           for(i=1;i<=ntr;i++){ 
              for(j=1;j<=ns;j++){ 

                 dummy_1[i][j] = (integrated_section_envelope[i][j] - integrated_sectiondata_envelope[i][j]) / (integrated_section_envelope[i][j]+EPS_LNORM6) * integrated_section_hilbert[i][j]; 

              } 
           } 
	
           calc_hilbert(dummy_1,dummy_2,ns,ntr); 

        }

        /* Log L2 envelope objective function*/
        if(ENV==2){

           for(i=1;i<=ntr;i++){
              for(j=1;j<=ns;j++){

                 dummy_1[i][j] = (log(integrated_sectiondata_envelope[i][j]+EPS_LNORM6) - log(integrated_section_envelope[i][j]+EPS_LNORM6)) / (integrated_section_envelope[i][j]*integrated_section_envelope[i][j]+EPS_LNORM6) * integrated_section_hilbert[i][j];

              }
           }

           calc_hilbert(dummy_1,dummy_2,ns,ntr);

        }
	

} /* end of if LNORM==6 */

/* Laplace-Fourier-domain residuals */
/* -------------------------------- */
if(INVMAT==1){
   laplace_fourier_res(sectiondata,section,sectiondiff,ntr,ntr_glob,ns,ishot,nsrc_glob,iter,recpos,recpos_loc,srcpos);
}	

/* Time-domain residuals */
/* --------------------- */
if(INVMAT!=1){
                  
/* calculate weighted data residuals and reverse time direction */
for(i=1;i<=ntr;i++){	
	
    if((TRKILL==1)&&(kill_vector[i]==1))
    continue;

    if(OFFSET_MUTE){
      
      /* calculate source and receiver positions */
      xr = recpos[1][recpos_loc[3][i]]*DH;
      xs = srcpos[1][ishot];
      yr = recpos[2][recpos_loc[3][i]]*DH;
      ys = srcpos[2][ishot];

      /* calculate absolute offset */
      offset = sqrt(((xs-xr)*(xs-xr))+((ys-yr)*(ys-yr)));

      if((OFFSET_MUTE==1)&&(offset>=OFFSETC)){continue;} /* mute far-offset data*/
      if((OFFSET_MUTE==2)&&(offset<=OFFSETC)){continue;} /* mute near-offset data*/

    }
	
    invtime=ns;
    
    intseis = 0.0;
    
    if (LNORM==5){
    	
    	abs_data=0.0;
	abs_synthetics=0.0;
	data_mult_synthetics=0.0;
	
        if((GRAD_FORM==2)||(GRAD_FORM==3)||(GRAD_FORM==4)){
	   for(j=1;j<=ns;j++){
	      intseis_data=sectiondata[i][j];
	      intseis_synthetics=section[i][j];
	      abs_data+=intseis_data*intseis_data;
	      abs_synthetics+=intseis_synthetics*intseis_synthetics;
	      data_mult_synthetics+=intseis_synthetics*intseis_data;
	   }
        }

        if(GRAD_FORM==1){
           for(j=1;j<=ns;j++){
              intseis_data=integrated_sectiondata[i][j];  
              intseis_synthetics=integrated_section[i][j];  
              abs_data+=intseis_data*intseis_data;  
              abs_synthetics+=intseis_synthetics*intseis_synthetics;  
              data_mult_synthetics+=intseis_synthetics*intseis_data;  
           }
        }

	abs_data=sqrt(abs_data);
	abs_synthetics=sqrt(abs_synthetics);
	
     }
	
      for(j=1;j<=ns;j++){
                       
                        if(j==1){sectiondata[i][j]=section[i][j];}
                        
			
			/* calculate L1 residuals */
			if(LNORM==1){
			if(((sectiondata[i][j]-section[i][j]))>0){signL1=1.0;}
			if(((sectiondata[i][j]-section[i][j]))<0){signL1=-1.0;}
			if(((sectiondata[i][j]-section[i][j]))==0){signL1=0.0;}
			
			sectiondiff[i][invtime]=signL1;
			}
			
			/* calculate L2 residuals */
			if(LNORM==2){

                          if(GRAD_FORM==1){
			     if(RTM==0){intseis += DT*(section[i][j]-sectiondata[i][j]);}
			     if(RTM==1){intseis += DT*(sectiondata[i][j]);}
                          }

                          if((GRAD_FORM==2)||(GRAD_FORM==3)||(GRAD_FORM==4)){
                             if(RTM==0){intseis = (section[i][j]-sectiondata[i][j]);}
                             if(RTM==1){intseis = sectiondata[i][j];}
                          }

			  sectiondiff[i][invtime]=intseis;

			}

			/* calculate Cauchy residuals */ 
		        if(LNORM==3){
			  sectiondiff[i][invtime]=((sectiondata[i][j]-section[i][j]))/(1+(((sectiondata[i][j]-section[i][j]))*((sectiondata[i][j]-section[i][j])))); 
			}
			
			/* calculate sech residuals */
			if(LNORM==4){
			  sectiondiff[i][invtime]=tanh(sectiondata[i][j]-section[i][j]);
			}  

                        /* calculate global correlation norm residuals and misfit function */
			if(LNORM==5){

                          if((GRAD_FORM==2)||(GRAD_FORM==3)||(GRAD_FORM==4)){
			    intseis_data = sectiondata[i][j];
			    intseis_synthetics = section[i][j];
                          }

                          if(GRAD_FORM==1){
                            intseis_data = integrated_sectiondata[i][j];
                            intseis_synthetics = integrated_section[i][j];
                          }

                          tmp = abs_synthetics*abs_synthetics*abs_synthetics*abs_data;
                          tmp1 = abs_synthetics*abs_data;

                          sectiondiff[i][invtime]=((intseis_synthetics*data_mult_synthetics)/(tmp+EPS_LNORM)) - (intseis_data/(tmp1+EPS_LNORM));

			}

                        if((LNORM==5)&&(swstestshot==1)){

                          tmp = abs_data*abs_synthetics;
                          L2-=(intseis_data*intseis_synthetics)/(tmp+EPS_LNORM);

	                }
			
                        /* Envelope objective functions after Chi, Dong and Liu (2014) */
                        if(LNORM==6){ 

                          if(ENV==1){ /* L2 envelope objective function */
                            sectiondiff[i][invtime]=integrated_section[i][j]/(integrated_section_envelope[i][j]+EPS_LNORM6)*(integrated_section_envelope[i][j]-integrated_sectiondata_envelope[i][j])-dummy_2[i][j];
                          }

                          if(ENV==2){ /* Log L2 envelope objective function */
                            sectiondiff[i][invtime]=dummy_2[i][j] - ((integrated_section[i][j]/(integrated_section_envelope[i][j]*integrated_section_envelope[i][j]+EPS_LNORM6))*(log(integrated_sectiondata_envelope[i][j]+EPS_LNORM6)-log(integrated_section_envelope[i][j]+EPS_LNORM6)));
                          }

                        } 

						
			invtime--;                                      /* reverse time direction */
       } 
}

   if(LNORM==6){
      /* low-pass filter envelope data */
      timedomain_filt(sectiondiff,FC,ORDER,ntr,ns,1);
      timedomain_filt(integrated_section_envelope,FC,ORDER,ntr,ns,1);
      timedomain_filt(integrated_sectiondata_envelope,FC,ORDER,ntr,ns,1);
   }

   if(TIMEWIN==3){
     stalta(sectiondiff, ntr, ns, picked_times, ishot);
     time_window(sectiondiff, picked_times, iter, ntr_glob,recpos_loc, ntr, ns, ishot);
   }

   if(TIMEWIN==4){
     time_window(sectiondiff, picked_times, iter, ntr_glob,recpos_loc, ntr, ns, ishot);
   }

   if((TIMEWIN==1)||(TIMEWIN==2)){  
     time_window(sectiondiff, picked_times, iter, ntr_glob,recpos_loc, ntr, ns, ishot);
   }    


} /* end of time domain residuals */

/* calculate objective functions */
/* ----------------------------- */

for(i=1;i<=ntr;i++){

      if((TRKILL==1)&&(kill_vector[i]==1))
      continue;

      if(OFFSET_MUTE){
      
        /* calculate source and receiver positions */
        xr = recpos[1][recpos_loc[3][i]]*DH;
        xs = srcpos[1][ishot];
        yr = recpos[2][recpos_loc[3][i]]*DH;
        ys = srcpos[2][ishot];

        /* calculate absolute offset */
       offset = sqrt(((xs-xr)*(xs-xr))+((ys-yr)*(ys-yr)));

       if((OFFSET_MUTE==1)&&(offset>=OFFSETC)){continue;} /* mute far-offset data*/
       if((OFFSET_MUTE==2)&&(offset<=OFFSETC)){continue;} /* mute near-offset data*/

     }	      

      invtime=ns;
      for(j=1;j<=ns;j++){


	 if((LNORM==2)&&(swstestshot==1)){
	   L2+=0.5 * sectiondiff[i][invtime]*sectiondiff[i][invtime]; 
	 }
	 
         /* Envelope objective functions after Chi, Dong and Liu (2014) */
	 if((LNORM==6) && (swstestshot==1)){ 

           if(ENV==1){ /* L2 envelope objective function */
              L2+=0.5 * (integrated_section_envelope[i][invtime]-integrated_sectiondata_envelope[i][invtime]) * (integrated_section_envelope[i][invtime]-integrated_sectiondata_envelope[i][invtime]); 
           }

           if(ENV==2){ /* Log L2 envelope objective function */
	      L2+=0.5 * (log(integrated_section_envelope[i][invtime]+EPS_LNORM6)-log(integrated_sectiondata_envelope[i][invtime]+EPS_LNORM6)) * (log(integrated_section_envelope[i][invtime]+EPS_LNORM6)-log(integrated_sectiondata_envelope[i][invtime]+EPS_LNORM6));
	   }


         }    
       
         invtime--;    /*reverse time direction */

	}
}

l2=L2;


if(TIMEWIN) free_vector(picked_times,1,ntr);

if((LNORM==5)&&(GRAD_FORM==1)){
  free_matrix(integrated_section,1,ntr,1,ns);
  free_matrix(integrated_sectiondata,1,ntr,1,ns);
}

if((LNORM==6)&&(GRAD_FORM==1)){
  free_matrix(integrated_section,1,ntr,1,ns);
  free_matrix(integrated_sectiondata,1,ntr,1,ns);
}

if(TRKILL){
  free_imatrix(kill_tmp,1,nsrc_glob,1,ntr_glob);
  free_ivector(kill_vector,1,ntr);
}

if(LNORM==6){ 
  free_matrix(integrated_sectiondata_envelope,1,ntr,1,ns); 
  free_matrix(integrated_section_envelope,1,ntr,1,ns); 
  free_matrix(integrated_section_hilbert,1,ntr,1,ns); 
  free_matrix(dummy_1,1,ntr,1,ns); 
  free_matrix(dummy_2,1,ntr,1,ns); 
} 
 
return l2;
} /* end of function */
