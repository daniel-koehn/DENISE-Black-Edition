/*------------------------------------------------------------------------
 *  Output of seismograms (PSV problem)
 * 
 *  Daniel Koehn
 *  Kiel, 25/04/2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void outseis_PSVfor(struct seisPSV *seisPSV, int *recswitch, int  **recpos, int  **recpos_loc, int ntr_glob, float ** srcpos, int ishot, int ns, int iter, FILE *FP){ 
		
        /* global variables */
	extern int SEISMO, MYID_SHOT, MODE;
	extern MPI_Comm SHOT_COMM;
	
        /* local variables */

	if ((SEISMO==1)&&(MODE==0)){
		
		catseis((*seisPSV).sectionvx, (*seisPSV).fulldata_vx, recswitch, ntr_glob, SHOT_COMM);
		catseis((*seisPSV).sectionvy, (*seisPSV).fulldata_vy, recswitch, ntr_glob, SHOT_COMM);
	
		if (MYID_SHOT==0){
		   saveseis_glob(FP,(*seisPSV).fulldata_vx,(*seisPSV).fulldata_vy,(*seisPSV).sectionp,(*seisPSV).sectioncurl,(*seisPSV).sectiondiv,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter);  
		}
	}

	if ((SEISMO==2)&&(MODE==0)){  

		catseis((*seisPSV).sectionp, (*seisPSV).fulldata_p, recswitch, ntr_glob, SHOT_COMM);

		if (MYID_SHOT==0){
		   saveseis_glob(FP,(*seisPSV).sectionvx,(*seisPSV).sectionvy,(*seisPSV).fulldata_p,(*seisPSV).sectioncurl,(*seisPSV).sectiondiv,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter);
		}
	}

	if ((SEISMO==3)&&(MODE==0)){  

		catseis((*seisPSV).sectioncurl, (*seisPSV).fulldata_curl, recswitch, ntr_glob, SHOT_COMM);
		catseis((*seisPSV).sectiondiv, (*seisPSV).fulldata_div, recswitch, ntr_glob, SHOT_COMM);

		if (MYID_SHOT==0){
		   saveseis_glob(FP,(*seisPSV).sectionvx,(*seisPSV).sectionvy,(*seisPSV).sectionp,(*seisPSV).fulldata_curl,(*seisPSV).fulldata_div,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter);
		}

	}


	if ((SEISMO==4)&&(MODE==0)){  

		catseis((*seisPSV).sectionvx, (*seisPSV).fulldata_vx, recswitch, ntr_glob, SHOT_COMM);
		catseis((*seisPSV).sectionvy, (*seisPSV).fulldata_vy, recswitch, ntr_glob, SHOT_COMM);
		catseis((*seisPSV).sectionp, (*seisPSV).fulldata_p, recswitch, ntr_glob, SHOT_COMM);
		catseis((*seisPSV).sectioncurl, (*seisPSV).fulldata_curl, recswitch, ntr_glob, SHOT_COMM);
		catseis((*seisPSV).sectiondiv, (*seisPSV).fulldata_div, recswitch, ntr_glob, SHOT_COMM);

		if (MYID_SHOT==0){
		   saveseis_glob(FP,(*seisPSV).fulldata_vx,(*seisPSV).fulldata_vy,(*seisPSV).fulldata_p,(*seisPSV).fulldata_curl,(*seisPSV).fulldata_div,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter);
		}         
	 
	}		
		
}
