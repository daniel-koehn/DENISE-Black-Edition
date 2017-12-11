/*------------------------------------------------------------------------
 *   write SH-seismograms to files 
 *   last update 09/12/2017, D. Koehn
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void saveseis_glob_SH(FILE *fp, float **sectionvz, int  **recpos, int  **recpos_loc, 
int ntr, float ** srcpos, int ishot, int ns, int iter){ 
		
	extern int SEISMO, SEIS_FORMAT, MYID, RUN_MULTIPLE_SHOTS, MODE;	
	extern char  SEIS_FILE_VX[STRING_SIZE]; 

        char vzf[STRING_SIZE];
        int nsrc=1;			

	if(MODE>0){
	  sprintf(vzf,"%s.shot%d.it%d",SEIS_FILE_VX,ishot,iter);
	}
		
	if(MODE==0){
	  sprintf(vzf,"%s.shot%d",SEIS_FILE_VX,ishot);
	}
	
	
	switch (SEISMO){
	case 1 : /* particle velocity (z-component) only */
		fprintf(fp," PE %d is writing %d seismograms (vz) to\n\t %s \n",MYID,ntr,vzf);
		outseis_glob(fp,fopen(vzf,"w"),1,sectionvz,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		break;			
        }     
}
