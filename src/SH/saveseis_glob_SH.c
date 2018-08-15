/*------------------------------------------------------------------------
 *   write SH-seismograms to files 
 *   last update 09/12/2017, D. Koehn
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void saveseis_glob_SH(FILE *fp, float **sectionvz, int  **recpos, int  **recpos_loc, 
int ntr, float ** srcpos, int ishot, int ns, int iter){ 
		
	extern int SEISMO, SEIS_FORMAT, MYID, RUN_MULTIPLE_SHOTS, MODE;	
	extern char  SEIS_FILE_VY[STRING_SIZE]; 

        char vyf[STRING_SIZE];
        int nsrc=1;			

	if(MODE>0){
	  sprintf(vyf,"%s.shot%d.it%d",SEIS_FILE_VY,ishot,iter);
	}
		
	if(MODE==0){
	  sprintf(vyf,"%s.shot%d",SEIS_FILE_VY,ishot);
	}
	
	
	switch (SEISMO){
	case 1 : /* particle velocity (y-component) only */
		fprintf(fp," PE %d is writing %d seismograms (vy) to\n\t %s \n",MYID,ntr,vyf);
		outseis_glob(fp,fopen(vyf,"w"),1,sectionvz,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		break;			
        }     
}
