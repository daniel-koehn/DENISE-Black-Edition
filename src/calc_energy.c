/*------------------------------------------------------------------------
 *   Calculate energy of measured data                                  
 *   last update 22/03/11, L. Rehor
 *  ----------------------------------------------------------------------*/
#include "fd.h"

double calc_energy(float **sectiondata, int ntr, int ns, float energy, int ntr_glob, int **recpos_loc, int nsrc_glob, int ishot){

/* declaration of variables */
extern float DT;
extern int MYID;
extern int TRKILL;
extern char TRKILL_FILE[STRING_SIZE];
int i, j;
float energy_dummy, intseis;
int umax=0, h;

/*extern int MYID;*/

/*printf("MYID: %d; ns=%d\n",MYID,ns);
printf("MYID: %d; ntr=%d\n",MYID,ntr);*/

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


for(i=1;i<=ntr;i++){
	if((TRKILL==1)&&(kill_vector[i]==1))
    	continue;	

        intseis = 0.0;
	for(j=1;j<=ns;j++){
	        intseis += sectiondata[i][j];
		energy += intseis*intseis*DT*DT;
	}
}

/*printf("MYID: %d; energy: %e\n",MYID,energy);*/

energy_dummy=energy;

return energy_dummy;

/* free memory for trace killing */
if(TRKILL){
free_imatrix(kill_tmp,1,ntr_glob,1,nsrc_glob);
free_ivector(kill_vector,1,ntr);
}
}
