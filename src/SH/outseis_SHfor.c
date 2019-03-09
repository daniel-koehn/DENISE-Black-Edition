/*------------------------------------------------------------------------
 *  Output of seismograms (SH problem)
 * 
 *  Daniel Koehn
 *  Kiel, 09/12/2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void outseis_SHfor(struct seisSH *seisSH, int *recswitch, int  **recpos, int  **recpos_loc, int ntr_glob, float ** srcpos, int ishot, int ns, int iter, FILE *FP){ 
		
        /* global variables */
	extern int SEISMO, MYID, MODE;	
	extern MPI_Comm SHOT_COMM;

        /* local variables */

	if ((SEISMO==1)&&(MODE==0)){
		
		catseis((*seisSH).sectionvz, (*seisSH).fulldata_vz, recswitch, ntr_glob, SHOT_COMM);
	
		if (MYID==0){
		   saveseis_glob_SH(FP,(*seisSH).fulldata_vz,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter);
		}
	}		
		
}
