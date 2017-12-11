/*
 * Allocate memory for seismograms of full data (SH problem) 
 *
 * Daniel Koehn
 * Kiel, 03/12/2017
 */

#include "fd.h"

void alloc_seisSHfull(struct seisSH *seisSH, int ntr_glob){

        /* global variables */
	extern int SEISMO, NT;

	/* local variables */	

	if(SEISMO==1){
  	   (*seisSH).fulldata_vz = matrix(1,ntr_glob,1,NT);
	}

}



