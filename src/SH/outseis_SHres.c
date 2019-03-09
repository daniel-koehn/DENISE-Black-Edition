/*------------------------------------------------------------------------
 *  Output of time reversed residual seismograms (SH problem)
 * 
 *  Daniel Koehn
 *  Kiel, 13/12/2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void outseis_SHres(struct seisSH *seisSH, struct seisSHfwi *seisSHfwi, int *recswitch, int  **recpos, int  **recpos_loc, int ntr_glob, float ** srcpos, int ishot, int ns, int nstage, FILE *FP){ 
		
        /* global variables */
	extern int QUELLTYPB, MYID;	
	
        /* local variables */

        if(QUELLTYPB){
   
           catseis((*seisSHfwi).sectionvzdiff, (*seisSH).fulldata_vz, recswitch, ntr_glob, MPI_COMM_WORLD);
      
           if (MYID==0){
	   	saveseis_glob_SH(FP,(*seisSH).fulldata_vz,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,nstage); 
           }
      
        }   		
}
