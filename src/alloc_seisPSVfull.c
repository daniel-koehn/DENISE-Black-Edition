/*
 * Allocate memory for seismograms of full data (PSV problem) 
 *
 * Daniel Koehn
 * Kiel, 24/04/2016
 */

#include "fd.h"

void alloc_seisPSVfull(struct seisPSV *seisPSV, int ntr_glob){

        /* global variables */
	extern int SEISMO, NT;

	/* local variables */	

	if(SEISMO){
  	   (*seisPSV).fulldata = matrix(1,ntr_glob,1,NT);
	}

	if(SEISMO==1){
  	   (*seisPSV).fulldata_vx = matrix(1,ntr_glob,1,NT);
  	   (*seisPSV).fulldata_vy = matrix(1,ntr_glob,1,NT);
	}

	if(SEISMO==2){
  	   (*seisPSV).fulldata_p = matrix(1,ntr_glob,1,NT);
	}

	if(SEISMO==3){
  	   (*seisPSV).fulldata_curl = matrix(1,ntr_glob,1,NT); 
  	   (*seisPSV).fulldata_div = matrix(1,ntr_glob,1,NT);
	}

	if(SEISMO==4){
  	   (*seisPSV).fulldata_vx = matrix(1,ntr_glob,1,NT);
  	   (*seisPSV).fulldata_vy = matrix(1,ntr_glob,1,NT);
  	   (*seisPSV).fulldata_p = matrix(1,ntr_glob,1,NT); 
  	   (*seisPSV).fulldata_curl = matrix(1,ntr_glob,1,NT);
  	   (*seisPSV).fulldata_div = matrix(1,ntr_glob,1,NT); 
	}

}



