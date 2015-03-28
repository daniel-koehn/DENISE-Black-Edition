/*
 *------------------------------------------------------------------------
 *
 *   modified energy ratio (MER) first arrival time picker
 *
 *------------------------------------------------------------------------
*/


#include "fd.h"

/* subroutine prototypes */
float mean(float *vec, int istart, int iend);
float maximumpos(float *a, int n);


void mer(float **sectionp, int ntr, int nst, float *picked_times, int ishot)
{
	
	extern float	DT;
	extern int POS[3];
	extern char PICKS_FILE[STRING_SIZE];

	int	ne;		/* sta time window in samples			*/
	int	ltawin;		/* lta time window in samples			*/
	float	*er;		/* energy ratio                                 */
	float   *er3;           /* energy ratio^3                               */
	int	maxpos=1;	/* position of maximum in stalta		*/
	float	*stmp, *ltmp;	/* temporary variables			        */
	int	j,k,i;
	float	*data, maxdata, nom, denom;
	char    jac[225];
        FILE    *fp1;

	/* define energy window */
        ne    = 40;	 /*40;*/
	
	/* allocate memory for vector variables */
	er = vector(1,nst);
	er3 = vector(1,nst);
	data = vector(1,nst);

        sprintf(jac,"%s_shot%d.%d%d",PICKS_FILE,ishot,POS[1],POS[2]);
        fp1=fopen(jac,"w");


	/* Main loop over traces */
	for(i=1;i<=ntr;i++) {
		/*****************/
		/* extract trace */
		/*****************/
		for (j=1;j<=nst;j++) {
			if (j==1) {
				data[1] = 0.0;
				continue;
			}
			data[j] = sectionp[i][j];
			if (j==2)  maxdata = data[j];
			else { if (maxdata<data[j]) maxdata = data[j]; }
		}
		
				
		/***********/
		/* MER */
		/***********/
		/* square of trace data */
		/*for (j=1;j<=nst;j++)  data[j] = data[j]*data[j]/(maxdata*maxdata);*/
		for (j=1;j<=nst;j++){data[j] = data[j]*data[j];}
		
		/* calculate nominator of energy ratio */
		for (j=1;j<=nst;j++){
			
			nom = 0.0;
			denom = 0.0;
			
			/* sum data over energy window */
			
			/* nominator */
			for (k=j-ne;k<=j;k++){
			  
			   if(k<=0){nom+=(data[1]+data[2])/2.0;}
			   if((k>0)&&(k<(nst+1))){nom+=data[k];}
			   if(k>=(nst+1)){nom+=(data[nst-1]+data[nst])/2.0;}
			   
			}
			 
			/* denominator */ 
			for (k=j;k<=j+ne;k++){ 
			                                                                            
			   if(k<=0){denom+=(data[1]+data[2])/2.0;} 
			   if((k>0)&&(k<(nst+1))){denom+=data[k];} 
			   if(k>=(nst+1)){denom+=(data[nst-1]+data[nst])/2.0;} 
			                                                                                                                                                                                       
			}
			
			if((nom==0)||(denom==0)){
			  er[j]=-1;
			}
			else{
			er[j] = nom/denom; /* energy ratio */
			}
			
			er3[j] = pow((fabs(data[j]*er[j])),3);
					  
		}
		
		maxpos = maximumpos(er3,nst);
		picked_times[i] = (float)maxpos*DT;
		
		/* output of picked times to file */
		fprintf(fp1,"%d \t %e \n",i,picked_times[i]);
	}
	
	fclose(fp1);
	
	free_vector(er,1,nst);
	free_vector(er3,1,nst);
	free_vector(data,1,nst);

}


/*========================================================================*/
/*                              SUBROUTINES                               */
/*========================================================================*/

/*float mean(float *vec, int istart, int iend)
{
	int i;
	float m=0;
	
	for (i=istart;i<=iend;i++) m += vec[i];
	m /= (iend-istart+1);
	
	return m;
}


float maximumpos(float *a, int n){
	float maxi=0.0;
	int j;
	int maxpos = 1;

	for (j=1;j<=n;j++)
	{	
		if (a[j]<0) continue;
		if (fabs(a[j])>maxi)
		{
			maxi = fabs(a[j]);
			maxpos = j;
		}
	}
		
	return maxpos;
}*/
