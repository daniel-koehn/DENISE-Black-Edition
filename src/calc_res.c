/*------------------------------------------------------------------------
 *   Calculate Data Residuals                                  
 *   last update 29/03/08, D.Koehn
 *  ----------------------------------------------------------------------*/
#include "fd.h"

double calc_res(float **sectiondata, float **section, float **sectiondiff, float **sectiondiffold, int ntr, int ns, int LNORM, float L2, int itest, int sws, int swstestshot, int ntr_glob, int **recpos, int **recpos_loc, float **srcpos, int nsrc_glob, int ishot, int iter){

/* declaration of variables */
extern float DT, DH, OFFSETC, FC, FC_START;
extern int REC1, REC2, MYID, ORDER;
extern int TRKILL, TIME_FILT, GRAD_FORM;
extern char TRKILL_FILE[STRING_SIZE];
extern int NORMALIZE, TIMEWIN, RTM, OFFSET_MUTE;
float RMS, RMS_obs, signL1, intseis;
int Lcount,i,j,invtime,k,h, umax=0;
float l2;
float abs_data, abs_synthetics, data_mult_synthetics, intseis_data, intseis_synthetics;
float intseis_section, intseis_sectiondata, offset, xr, yr, xs, ys;
float *picked_times=NULL, eps;
float **integrated_section=NULL, **integrated_sectiondata=NULL;
float **sectiondata_envelope=NULL, **section_envelope=NULL, **hilbert=NULL;
float **dummy_1=NULL, **dummy_2=NULL; 
float EPS_LNORM;

EPS_LNORM=1e-3;
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

if(LNORM==6){

	integrated_section = matrix(1,ntr,1,ns);
	integrated_sectiondata = matrix(1,ntr,1,ns);

/* Integration before TIMEWIN of measured and synthetic data  */
intseis_section = 0.0;
for(i=1;i<=ntr;i++){
      	for(j=1;j<=ns;j++){
		intseis_section += section[i][j];
		integrated_section[i][j]=intseis_section*DT;
		}
	}
intseis_sectiondata = 0.0;
for(i=1;i<=ntr;i++){
      	for(j=1;j<=ns;j++){
		intseis_sectiondata += sectiondata[i][j];
		integrated_sectiondata[i][j]=intseis_sectiondata*DT;
		}
	}
} /* end of if LNORM==6 */

if((TIMEWIN==1)||(TIMEWIN==2)){
  
  if(LNORM==6){
    time_window(integrated_sectiondata, picked_times, iter, ntr_glob,recpos_loc, ntr, ns, ishot);
    time_window(integrated_section, picked_times, iter, ntr_glob,recpos_loc, ntr, ns, ishot);
  }
  
  time_window(sectiondata, picked_times, iter, ntr_glob,recpos_loc, ntr, ns, ishot);
  time_window(section, picked_times, iter, ntr_glob,recpos_loc, ntr, ns, ishot);
}
                  
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
	
	for(j=1;j<=ns;j++){
		intseis_data=sectiondata[i][j];
		intseis_synthetics=section[i][j];
		abs_data+=intseis_data*intseis_data;
		abs_synthetics+=intseis_synthetics*intseis_synthetics;
		data_mult_synthetics+=intseis_synthetics*intseis_data;
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

                          if((GRAD_FORM==2)||(GRAD_FORM==3)){
                             if(RTM==0){intseis = section[i][j]-sectiondata[i][j];}
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
			  intseis_data = sectiondata[i][j];
			  intseis_synthetics = section[i][j];

                          if((fabs(abs_synthetics)>eps)&&(fabs(abs_data)>eps)){
			    sectiondiff[i][invtime]=((intseis_synthetics*data_mult_synthetics)/(abs_synthetics*abs_synthetics*abs_synthetics*abs_data)) - (intseis_data/(abs_synthetics*abs_data));
                          }else{sectiondiff[i][invtime]=0.0;}

			}

                        if((LNORM==5)&&(swstestshot==1)){
                          if((fabs(abs_synthetics)>eps)&&(fabs(abs_data)>eps)){
	                      L2-=(intseis_data*intseis_synthetics)/(abs_data*abs_synthetics);
                          }
	                }
				
			if(LNORM==6){
			sectiondiff[i][invtime]=integrated_section[i][j]-integrated_sectiondata[i][j];
			}
						
			invtime--;                                      /* reverse time direction */
       } 
}

if(TIMEWIN==3){
  stalta(sectiondiff, ntr, ns, picked_times, ishot);
  time_window(sectiondiff, picked_times, iter, ntr_glob,recpos_loc, ntr, ns, ishot);
}

if(TIMEWIN==4){
  time_window(sectiondiff, picked_times, iter, ntr_glob,recpos_loc, ntr, ns, ishot);
}    


/* calculate misfit functions */
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
	   L2+=sectiondiff[i][invtime]*sectiondiff[i][invtime]; 
	 }
					
	 if((LNORM==6)&&(swstestshot==1)){
	   L2+=sectiondiff[i][invtime]*sectiondiff[i][invtime];
	 }
         
         invtime--;    /*reverse time direction */

	}
}

l2=L2;


if(TIMEWIN) free_vector(picked_times,1,ntr);
if(LNORM==6){
free_matrix(integrated_section,1,ntr,1,ns);
free_matrix(integrated_sectiondata,1,ntr,1,ns);
}

if(TRKILL){
free_imatrix(kill_tmp,1,nsrc_glob,1,ntr_glob);
free_ivector(kill_vector,1,ntr);
} 

return l2;
} /* end of function */
