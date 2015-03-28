/*------------------------------------------------------------------------
 *   write seismograms to files 
 *   last update 19/01/02, T. Bohlen
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void saveseis_test(FILE *fp, float **sectionp,
int  **recpos, int  **recpos_loc, 
int ntr, float ** srcpos, int ishot,int ns){ 
		
	extern int SEISMO, SEIS_FORMAT, MYID;	
	extern char  SEIS_FILE_P[STRING_SIZE];

        char pf[STRING_SIZE];
        int nsrc=1;		
	
        sprintf(pf,"su/test/TEST.shot%d.%d",ishot,MYID);

	outseis(fp,fopen(pf,"w"), 0, sectionp,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT);

}
