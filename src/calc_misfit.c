/*------------------------------------------------------------------------
 *   Calculate Misfit                                  
 *   last update 18/04/11, L. Rehor
 *  ----------------------------------------------------------------------*/
#include "fd.h"

double calc_misfit(float **sectiondiff, int ntr, int ns, int LNORM, float L2, int ntr_glob, int **recpos_loc, int nsrc_glob, int ishot){

/* declaration of variables */
extern float DT;
extern int MYID;
extern int TRKILL;
extern char TRKILL_FILE[STRING_SIZE];
int i,j;
float l2;
float L2_dummy;
int umax=0, h;
	
/* declaration of variables for trace killing */
int ** kill_tmp, *kill_vector;
char trace_kill_file[STRING_SIZE];	
FILE *ftracekill;

if(TRKILL){
	kill_tmp = imatrix(1,ntr_glob,1,nsrc_glob);
	kill_vector = ivector(1,ntr);

	ftracekill=fopen(TRKILL_FILE,"r");

	if (ftracekill==NULL) err(" Trace kill file could not be opened!");

	for(i=1;i<=ntr_glob;i++){
		for(j=1;j<=nsrc_glob;j++){
			fscanf(ftracekill,"%d",&kill_tmp[i][j]);
		}
	}

	fclose(ftracekill);

	h=1;
	for(i=1;i<=ntr;i++){
	   kill_vector[h] = kill_tmp[recpos_loc[3][i]][ishot];
	   h++;
	}
} /* end if(TRKILL)*/

/* calculate misfit for misfit function 1-6 */
for(i=1;i<=ntr;i++){
      if((TRKILL==1)&&(kill_vector[i]==1))
      continue;	
 
    
      for(j=1;j<=ns;j++){
 	                        
	/* calculate norm */
	if((LNORM==2) || (LNORM==6)){
	   L2+=sectiondiff[i][j]*sectiondiff[i][j]; 
	}
			
      } /* end of loop over time samples */
      

} /* end of loop over traces */


l2=L2;
return l2;

/* free memory for trace killing */
if(TRKILL){
free_imatrix(kill_tmp,1,ntr_glob,1,nsrc_glob);
free_ivector(kill_vector,1,ntr);
}

}
