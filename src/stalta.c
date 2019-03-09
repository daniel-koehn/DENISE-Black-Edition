/*
 *------------------------------------------------------------------------
 *
 *   STA/LTA first arrival time picker
 *
 *------------------------------------------------------------------------
*/


#include "fd.h"

/* subroutine prototypes */
void cumsum(float *vec, float *cumsum2, int n);
float maximumpos(float *a, int n);

void stalta(float **sectionp, int ntr, int nst, float *picked_times, int ishot)
{
	
	extern float DT;
	extern int POS[3], MYID;
	extern char PICKS_FILE[STRING_SIZE2];

	int	stawin;		/* sta time window in samples			*/
	int	ltawin;		/* lta time window in samples			*/
	float	*sta;		/* average of samples over window stawin	*/
	float	*lta;		/* average of samples over window ltawin	*/
	float	*stalta;	/* sta/lta					*/
        float   TR;             /* Value of STA/LTA ratio that triggers         */
        float   DTR;            /* Value of STA/LTA ratio that untriggers       */
	int	staltawin;	/* search max(sta/lta) in (1,...,staltawin)	*/
	int	maxpos=1;	/* position of maximum in stalta		*/
	float	*stmp, *ltmp;	/* temporary variables			*/
	int	j,k,i,ntshift,h;
	float	*data, *data1, *trend, **A, **A1, **A1inv, *coeff;
        float   *cumsuml, *cumsums;
        float   maxdata, maxstalta, tshift, tshift1;
        float   sumtmp;
	char    jac[225];
        FILE *fp1;

	/* define time windows */
        stawin    = ceil(0.007/DT);
	ltawin    = ceil(0.08/DT); 
	staltawin = nst;
        TR = 3.0;
        DTR = 1.0;
	tshift = 0.0;  /* 0.2 */  /* timeshift [s] to avoid picking errors at small offsets */
	tshift1 = 0.0; /* 0.2 */  /* arbitarily timeshift [s] for picked arrival times */
	ntshift = floor(tshift/DT);  /* timeshift [timesamples] */
	
	/* allocate memory for vector variables */
	lta = vector(1,nst);
	sta = vector(1,nst);
	ltmp = vector(1,nst-ltawin+1);
	stmp = vector(1,nst-stawin+1);
        cumsuml = vector(1,nst-ltawin+1);
        cumsums = vector(1,nst-stawin+1);
	stalta = vector(1,nst);
	data = vector(1,nst);
        data1 = vector(1,nst);
        trend = vector(0,nst);
        coeff = vector(1,2);
        A = matrix(0,nst,0,2);
    
        sprintf(jac,"%s_shot_%d",PICKS_FILE,ishot);
        fp1=fopen(jac,"w");

	/* Main loop over traces */
	for(i=1;i<=ntr;i++) { 

		/*****************/
		/* extract trace */
		/*****************/
                h=1;
		for (j=1;j<=nst;j++) {
			data[h] = sectionp[i][j];
                h++;
		}
		
		/****************************/
		/* Apply time shift to trace*/
		/****************************/
		/*for (j=1;j<=nst;j++) {
		    
		    if((j+ntshift)<=nst){data1[j+ntshift]=data[j];}
		    
		}*/

                /****************************/
		/* Detrend the time series  */
                /****************************/
                h=1;
                for (j=0;j<=nst-1;j++) {
		    
                    trend[j] = data[h];	
                    A[j][0] = j;
	            A[j][1] = 1.0;
		    h++;    
		}

                /* solve matrix equation */
                householder(nst,2,A,trend);

                
		/*for (h=1;h<=2;h++){
                     for (j=1;j<=2;j++){
                     
                         printf("A1 = %e \t h=%d \t j=%d \n ",A1[h][j],h,j);
                     
                     }
                }*/   
                
                /*printf("coeff = %e \t %e \t %d \t %d \n",trend[0],trend[1],i,MYID);*/

                /* calculate and remove linear trend */
                for (j=1;j<=nst;j++) {
		    
                    data[j] = data[j] - trend[0]*j - trend[1];	
		    
		}
				
		/***********/
		/* STA/LTA */
		/***********/
		/* square of trace data */
		for (j=1;j<=nst;j++)  data[j] = sqrt(data[j]*data[j]);
		
		/* apply STA/LTA-algorithm to trace */
                
                /* long-term average*/
                sumtmp=0.0;
                for (j=1;j<=ltawin;j++){
                    sumtmp +=data[j];
                }

                ltmp[1] = sumtmp;
                h = ltawin+1;
                k = 1;
                for (j=2;j<=(nst-ltawin+1);j++){
                    ltmp[j] = data[h]-data[k];
                    h++;
                    k++;
                }

                cumsum(ltmp,cumsuml,nst-ltawin+1);

                /* short-term average*/
                sumtmp=0.0;
                for (j=1;j<=stawin;j++){
                    sumtmp +=data[j];
                }

                stmp[1] = sumtmp;
                /*if((MYID==0)&&(i==1)){
                printf("stawin= %d \t stmp[1]=%e \n",stawin,sumtmp);}*/
                h = stawin+1;
                k = 1;
                for (j=2;j<=(nst-stawin+1);j++){
                    stmp[j] = data[h]-data[k];
                    h++;
                    k++;
                }

                cumsum(stmp,cumsums,nst-stawin+1);

                /* calculate short-term / long-term average */
                h=1;
                k=1;
                for (j=1;j<=nst;j++){

                    if(j<=ltawin-1){lta[j]=cumsuml[1]/ltawin;}
                    if(j>=ltawin){
                      lta[j]=cumsuml[h]/ltawin;
                      h++;
                    }

                    if(j<=stawin-1){sta[j]=cumsums[1]/stawin;}
                    if(j>=stawin){
                      sta[j]=cumsums[k]/stawin;
                      k++;
                    }

                    stalta[j] = sta[j]/lta[j];

                }
                
  
		/*for (j=1;j<=nst;j++) {
			if (lta[j] == 0 || sta[j] == 0)
				stalta[j] = -1; *//*setnan(stalta[j]);*/
		/*	else
				stalta[j] = sta[j]/lta[j];	

		}*/
		maxpos = maximumpos(stalta,nst);
                /*maxpos = maxpos - (int)(((float)(ltawin-stawin))/2);*/
		/*picked_times[i] = (float)maxpos*DT-tshift+tshift1;*/
   		picked_times[i] = (float)maxpos*DT;


		/* output of picked times to file */
	        fprintf(fp1,"%e \t %e \n",picked_times[i],picked_times[i]);
	
              


	}	

        fclose(fp1);
	
	free_vector(lta,1,nst);
	free_vector(sta,1,nst);
        free_vector(cumsuml,1,nst-ltawin+1);
        free_vector(cumsums,1,nst-stawin+1);
	free_vector(ltmp,1,nst-ltawin+1);
	free_vector(stmp,1,nst-stawin+1);
	free_vector(stalta,1,nst);
	free_vector(data,1,nst);
	free_vector(trend,0,nst);
	free_vector(coeff,1,2);
        free_matrix(A,0,nst,0,2);

}


/*========================================================================*/
/*                              SUBROUTINES                               */
/*========================================================================*/

void cumsum(float *vec, float *cumsum2, int n)
{
	int i;
	
        cumsum2[1] = vec[1];
	for (i=2;i<=n;i++){
	    cumsum2[i] = cumsum2[i-1] + vec[i];
        }

}


float maximumpos(float *a, int n){
	float maxi=0.0;
	int j;
	int maxpos = 1;

	for (j=1;j<=n;j++)
	{
		/* if (isnan(a[j])) continue; */
		/* if (!finite(a[j])) continue; */
		if (a[j]<0) continue;
		if (fabs(a[j])>maxi)
		{
			maxi = fabs(a[j]);
			maxpos = j;
		}
	}
		
	return maxpos;
}
